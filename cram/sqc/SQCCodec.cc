/** @file SQCCodec.cc
 *  @brief This file contains the implementation of the SQCCodec class.
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <map>
#include <string.h>
#include <utility>

#include "CABAC_ArithmeticEncoder.h"
#include "CABAC_ArithmeticDecoder.h"
#include "CABAC_BitstreamFile.h"
#include "ContextModel.h"

#include "SQCCodec.h"
#include "htslib/hts_log.h"

namespace sqc {

// BLOCK_GROW and BLOCK_RESIZE operators are defined in cram_encode.c
#ifndef BLOCK_GROW
/* Block size and data pointer. */
#define BLOCK_SIZE(b) ((b)->byte)
/* Request block to be at least 'l' bytes long */
#define BLOCK_RESIZE(b,l)                                       \
		do {                                                        \
			while((b)->alloc <= (l)) {                              \
				(b)->alloc = (b)->alloc ? (b)->alloc*1.5 : 1024;    \
				(b)->data = (unsigned char*)realloc((b)->data, (b)->alloc);         \
			}                                                       \
		} while(0)
/* Ensure the block can hold at least another 'l' bytes */
#define BLOCK_GROW(b,l) BLOCK_RESIZE((b), BLOCK_SIZE((b)) + (l))
#endif

int SQCCodec::qualityValueOffset_ = 33;
static uint8_t bit_mask[8] = {0b1,0b10,0b100,0b1000,0b10000,0b100000,0b1000000,0b10000000}; 

int header_buffer_append(uint32_t val, int len, int &bit, uint8_t **ptr)
{
	while (len)
	{
		int seg = min(bit,len); 
		(**ptr) |= (uint8_t)(((val >> (len-seg)) & ((1<<seg) -1 )) << (8 - bit));
		len -= seg; 
		bit -= seg; 
		if (!bit){
			(*ptr)++;
			(**ptr) = 0; 
			bit = 8;
		}
	}
	return 0;
}
int header_buffer_fetch(int len, int &bit, uint8_t **ptr)
{
	int val = 0;
	int seg = 0; 
	while (len)
	{
		seg = min(bit,len); 
		val <<= seg; 
		val |= ((**ptr) >> (8 - bit)) & ((1 << seg) - 1);
		bit -= seg;
		len -= seg; 
		if (!bit){
			(*ptr)++;
			bit = 8;
		}
	}
	return val; 
}

SQCCodec::SQCCodec()
    : nrMappedRecords_(0),
      nrUnmappedRecords_(0),
	  samPileupDeque_(),
      samRecordDeque_()
{
	samPileupDeque_.clear(); 
	samRecordDeque_.clear();
    samRecordDeque_.reserve(SEQS_PER_SLICE); // as defined in cram_structs.h

    maxSegID_ = 0;
    for (int d = 0; d < 4; d++) {
    	for (int i = 0; i < MAX_SEGMENTS; i++) {
    		qualitySegMin_[d][i] = 63;
		}
    }

    compressedQualSize_ = 0;
    uncompressedQualSize_ = 0;
    headerSize_ = 0;

    refOffsetVec_.clear();
	genoUncVec_ = NULL;
	genoUncVecSize_ = 0;
}

SQCCodec::~SQCCodec(void) {}

void SQCCodec::updateQualitySegMinMax(const bam_seq_t *b, sqc_option opt)
{
	int32_t pos = 0;
	size_t seg_id = 0, seg_pos = 0; 
	uint16_t d = (bam_is_rev(b) ? 2 : 0) + ((b->core.flag & BAM_FREAD1) != 0);
	uint8_t *qual = bam_get_qual(b);
	uint8_t *base = bam_get_seq(b); 

	while (pos < b->core.l_qseq) {
		uint8_t q = qual[pos];
		uint8_t c = "=ACMGRSVTWYHKDBN"[bam_seqi(base, pos)];
		
		if ((q < qualitySegMin_[d][seg_id]) && (c != 'N'))
			qualitySegMin_[d][seg_id] = q;
		if (q > qualitySegMax_[d][seg_id])
			qualitySegMax_[d][seg_id] = q;

		if (qualitySegValueValid_[d][seg_id][q] == 0)     
			qualitySegValueValid_[d][seg_id][q] = 1;

		pos++; 
		seg_pos++;
		if (seg_pos >= opt.segsize) {
			seg_id++;
			seg_pos = 0;
		}
	}
	
	seg_id = ((b->core.l_qseq-1)/opt.segsize)+1;
	if (maxSegID_ < seg_id) maxSegID_ = seg_id; 
}

uint16_t SQCCodec::updateQualitySegBPMax(bool enc)
{
	uint16_t d, seg_id;
	uint16_t q, qmin, qmax, q1, max_bp = 0;

	for (d = 0; d < 4; d++) {
		for (seg_id = 0; seg_id < maxSegID_; seg_id++) {
			qmin = qualitySegMin_[d][seg_id];
			qmax = qualitySegMax_[d][seg_id]; 

			if (qmin > qmax) continue; // may happen when data does not cover all types of alignments, e.g., single-end

			q1 = 0;
			if (enc) { // mapping at encoder side
				for (q = qmin; q <= qmax; q++) {
					if (qualitySegValueValid_[d][seg_id][q]) {
						qualitySegValueMap_[d][seg_id][q] = q1;
						q1++;
					}
				}
			} else { // reverse mapping at decoder side
				for (q = qmin; q <= qmax; q++) {
					if (qualitySegValueValid_[d][seg_id][q]) {
						qualitySegValueMap_[d][seg_id][q1] = q;
						q1++;
					}
				}
			}
			q1--; 
			qualitySegMax_c[d][seg_id] = q1;   // max quality after range compression
			qualitySegBPMax_[d][seg_id] = maxBitplane(0, q1);

			if (max_bp < qualitySegBPMax_[d][seg_id])
				max_bp = qualitySegBPMax_[d][seg_id];
		}
	}
	// std::cout << max_bp << std::endl;
	return max_bp;
}

void SQCCodec::initRefOffsetVec(int nref)
{
	refOffsetVec_.resize(nref, -1);
}

void SQCCodec::addUnmappedRecordToBlock(const SAMRecord &samRecord, sqc_option opt)
{
    samRecordDeque_.push_back(samRecord);
//    samRecord.printLong();

    uncompressedQualSize_ += samRecord.seq.length();
    nrUnmappedRecords_++;
}

void SQCCodec::addMappedRecordToBlock(SAMRecord &samRecord, sqc_option opt)
{
#if debugging_code
	if (refOffsetVec_.empty()) {
		hts_log_error("SQC: refOffsetVec_ is empty");
	}
#endif
	if (refOffsetVec_[samRecord.ref_id] < 0) // if this ref_id is not seen before
		refOffsetVec_[samRecord.ref_id] = samPileupDeque_.posMax(); 

	samRecord.pos += refOffsetVec_[samRecord.ref_id];
	samRecord.posMin += refOffsetVec_[samRecord.ref_id];
	samRecord.posMax += refOffsetVec_[samRecord.ref_id];

    if (nrMappedRecords() == 0) {
        samPileupDeque_.setPosMin(samRecord.posMin);
        samPileupDeque_.setPosMax(samRecord.posMax);
    }
    else if (samRecord.posMax > samPileupDeque_.posMax()) {
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    samRecord.addToPileupQueue(samPileupDeque_, samRecord.qual);
    samRecordDeque_.push_back(samRecord);
//    samRecord.printLong();

    uncompressedQualSize_ += samRecord.seq.length();
    nrMappedRecords_++;
}

void SQCCodec::encodeBlock(cram_block *b)
{
	encodeMappedQual();
	writeBlock(b);
}

void SQCCodec::decodeBlock(cram_block *b, cram_block *bout)
{
	if (!b) return;

	// std::cout << samPileupDeque_.posMin() << ", " << samPileupDeque_.posMax() << std::endl;
	uint8_t *p = b->data;
	uint32_t nSamRecords = samRecordDeque_.size();

	compressedQualSize_ = b->comp_size;
	headerSize_ = b->comp_size - compressedQualSize_;

#if debugging_code
	printf("%d, %d\n", b->comp_size, b->uncomp_size);
	printf("%03u, ", compressedQualSize_);
	for (uint32_t j = 0; j < (b->comp_size - headerSize_); j++)
		printf("%03u, ", (uint32_t)(p[j]));
	printf("\n");
#endif

	decodeMappedQual(p);

#if debugging_code
	for (uint32_t i = 0; i < nSamRecords; i++) {
		const SAMRecord &sr = samRecordDeque_.at(i);
		for (int j = 0; j < sr.qual.size(); j++)
			std::cout << (uint8_t)(sr.qual[j] + 33);
		std::cout << std::endl;
	}
#endif
	// replace quality strings with decoded values
	p = bout->data;
	for (uint32_t i = 0; i < nSamRecords; i++) {
		const SAMRecord &sr = samRecordDeque_[i];
		memcpy(p, sr.qual.c_str(), sr.readLen);
		p += sr.readLen;
	}
}

size_t SQCCodec::writeBlock(cram_block *b)
{
	uint8_t *header_buf = new uint8_t[DATA_BUF_SIZE];
	int header_len = writeQualitySegMinMaxTable(header_buf);
	
	b->comp_size = compressedQualSize_ + header_len; 
	b->method = AGI_SQC;
	b->uncomp_size = uncompressedQualSize_;
	
	if (b->comp_size > (int32_t)(b->alloc))
		BLOCK_GROW(b, b->comp_size);

	// write header bytes
	memcpy(b->data, header_buf, header_len);
	// write data bytes
	memcpy(b->data+header_len, buf_.dat, compressedQualSize_);

#if debugging_code
	printf("%d, %d\n", b->comp_size, b->uncomp_size);
	printf("%03u, ", compressedQualSize_);
	p = b->data;
	for (uint32_t j = 0; j < b->comp_size; j++)
		printf("%03u, ", (uint32_t)(p[j]));
	printf("\n");
#endif
	delete[] header_buf;
	return b->comp_size; // return the total number of bytes written
}

int16_t SQCCodec::writeQualitySegMinMaxTable(uint8_t * buf)
{
	uint32_t b,d, i, q;
	int bit=8; 
	uint8_t *buf0 = buf; 
	*buf = 0; 
	header_buffer_append(maxSegID_,6,bit,&buf0); 
	header_buffer_append(opt_.segsize, 8, bit, &buf0); 
	header_buffer_append(opt_.sqc_cbr_bp,3,bit,&buf0); 

#if debugging_code
	std::cout << "writeQualitySegMinMaxTable begins: " << maxSegID_ << ", " << opt_.segsize << std::endl;
#endif
	for (d = 0; d < 4; d++) {
		for (i = 0; i < maxSegID_; i++) {
			header_buffer_append(qualitySegMin_[d][i], 6, bit, &buf0); 
			header_buffer_append(qualitySegMax_[d][i], 6, bit, &buf0); 
			for (b=offset_table[opt_.sqc_cbr_bp-1];b<64;b++){
				uint32_t val = cbQualRecOffset_[d][i][b];
				header_buffer_append(val, 6, bit, &buf0); 
				// std::cout << "d:" << d << " seg:" << i << " b:" << b << " val:" << val  << std::endl; 
			}
			for (q = qualitySegMin_[d][i]; q <= qualitySegMax_[d][i]; q++)
				header_buffer_append(qualitySegValueValid_[d][i][q], 1, bit, &buf0); 
		}
	}

#if debugging_code
	for (d = 0; d < 4; d++) {
		for (i = 0; i < maxSegID_; i++) {
			std::cout << qualitySegMin_[d][i] << ", ";
		}                 
		std::cout << std::endl;
	}
	for (d = 0; d < 4; d++) {
		for (i = 0; i < maxSegID_; i++) {
			std::cout << qualitySegMax_[d][i] << ", ";
		}
		std::cout << std::endl;
	}
// 	for (d = 0; d < 4; d++) {
// 		for (i = 0; i < maxSegID_; i++) {
// 			for (q = qualitySegMin_[d][i]; q <= qualitySegMax_[d][i]; q++) {
// 				if (qualitySegValueValid_[d][i][q]) std::cout << "1";
// 				else std::cout << "0";
// 			}
// 			std::cout << std::endl;
// 		}
// 	}
	std::cout << "writeQualitySegMinMaxTable ends\n";
#endif
	int hd = buf0-buf+((8==bit)?0:1);
	headerSize_ += hd; 
	return hd;
}


int16_t SQCCodec::readQualitySegMinMaxTable(uint8_t *p)
{
	unsigned int b, d, i, q;
	int bit = 8; 
	uint8_t *p0 = p; 
	
	maxSegID_ = header_buffer_fetch(6,bit,&p0); 
	// std::cout << "maxSegID_:" << maxSegID_ ; 
	opt_.segsize = header_buffer_fetch(8,bit,&p0); 
	// std::cout << "segsize" << opt_.segsize;  

#if debugging_code
	std::cout << "readQualitySegMinMaxTable begins: " << maxSegID_ << ", " << opt_.segsize << std::endl;
#endif

	opt_.sqc_cbr_bp = header_buffer_fetch(3,bit,&p0);  // lowest bitplane that we use context based decoding offset  1...6
	// printf("in SQCCodec::readQualitySegMinMaxTable sqc_cbr_bp = %d\n", opt_.sqc_cbr_bp);
	for (d = 0; d < 4; d++) {
		for (i = 0; i < maxSegID_; i++) {
			qualitySegMin_[d][i] = header_buffer_fetch(6,bit,&p0); 
			qualitySegMax_[d][i] = header_buffer_fetch(6,bit,&p0); 
			for (b=offset_table[opt_.sqc_cbr_bp-1];b<64; b++) 
			{
				unsigned int val; 
				val =  header_buffer_fetch(6,bit,&p0);
				// std::cout << "d:" << d << " seg:" << i << " b:" << b << " val:" << val  << std::endl; 
				cbQualRecOffset_[d][i][b] = val*1.0; 
			}
			for (q = qualitySegMin_[d][i]; q <= qualitySegMax_[d][i]; q++)
				 qualitySegValueValid_[d][i][q] = header_buffer_fetch(1,bit,&p0);
		}
	}

#if debugging_code
	for (d = 0; d < 4; d++) {
		for (i = 0; i < maxSegID_; i++) {
			std::cout << qualitySegMin_[d][i] << ", ";
		}
		std::cout << std::endl;
	}
	for (d = 0; d < 4; d++) {
		for (i = 0; i < maxSegID_; i++) {
			std::cout << qualitySegMax_[d][i] << ", ";
		}
		std::cout << std::endl;
	}

// 	for (d = 0; d < 4; d++) {
// 		for (i = 0; i < maxSegID_; i++) {
// 			for (q = qualitySegMin_[d][i]; q <= qualitySegMax_[d][i]; q++) {
// 				if (qualitySegValueValid_[d][i][q]) std::cout << "1";
// 				else std::cout << "0";
// 			}
// 			std::cout << std::endl;
// 		}
// 	}
#endif
	return p0-p+((8==bit)?0:1);
}

void SQCCodec::truncateBitstream(cram_block *b, double bpq_out)
{
	int32_t nRemaining = b->uncomp_size; // number of quality values
	int32_t nMaxBytes = (bpq_out >= 0) ? (uint32_t)(floor(nRemaining * bpq_out / 8.0)) : 0xffffffff; // bytes per quality
	
	buf_.init();
	if (nMaxBytes > b->comp_size) {
		std::cout << "warning: bitstream not truncated" << std::endl;
	} else {
		buf_.setBuf(b->data, nMaxBytes);
	}
}

int SQCCodec::getTruncatedBitstream(uint8_t *b) const
{
	return buf_.getBuf(b);
}

int SQCCodec::setTruncatedBitstream(const uint8_t *b, int len)
{
	return buf_.setBuf(b, len);
}

void SQCCodec::encodeMappedQual()
{
	uint8_t bi;
	uint16_t bp_max, d, bp, max_bs, prev_bp, prev_q = 0; 
	uint32_t nRemaining, nMaxBytes;
	uint8_t *q; 
	
	nMaxBytes = (opt_.sqc_bpq >= 0) ? (uint32_t)(floor(uncompressedQualSize_ * opt_.sqc_bpq / 8.0)) : 0xffffffff; // bytes per quality
	// std::cout << "bpq: " << opt_.sqc_bpq << ", nMaxBytes: " << nMaxBytes << std::endl;
	buf_.init();
	buf_.init_ac(ctxt); 

	genoUncVec_ = (GenoUnc *)calloc(uncompressedQualSize_, sizeof(GenoUnc));
	if (!genoUncVec_)
		hts_log_error("SQC: memory failure");

	uint32_t l = buf_.getlen();

	bp_max = updateQualitySegBPMax(true);

	if (samRecordDeque_.empty()) {
		std::cout << "samRecordDeque_ is empty\n";
		return;
	}
	// calculate bit shifting for all quality values
	max_bs = encodeQualBitshifting(buf_) + bp_max;

	q = (uint8_t *)calloc(genoUncVecSize_, sizeof(uint8_t)); 
	if (!q)
		hts_log_error("SQC: memory failure");

	for (size_t i = 0; i < genoUncVecSize_; i++) {
		GenoUnc &gu = genoUncVec_[i];
		const SAMRecord &sr = samRecordDeque_[gu.read_idx];
		int q0 = sr.qual[gu.qual_idx]-SQCCodec::qualityValueOffset_; 
		q[i] = qualitySegValueMap_[sr.cidx][gu.seg_id][q0];
		for (int offset_bp = 0; offset_bp < 6; offset_bp++) /* build up offset table for partially decoded quality values */
		{
			size_t offset_ctx = (q[i]>>(offset_bp + 1))+SQCCodec::offset_table[offset_bp];  /* we start using context based offset when the quantization interval >= 2^3-1 */
			if (cbQualRecOffset_[sr.cidx][gu.seg_id][offset_ctx] < 1.0)
				cbQualRecOffset_[sr.cidx][gu.seg_id][offset_ctx] = q0;
			else
				cbQualRecOffset_[sr.cidx][gu.seg_id][offset_ctx] = 0.96*cbQualRecOffset_[sr.cidx][gu.seg_id][offset_ctx] + 0.04*q0;  /* mean quality value for each context */
		}
	}

	// ------------------------------------------------------------------------------------------
	// scan all the quality values and encode via bit shifting values
	nRemaining = uncompressedQualSize_;
	prev_bp = max_bs + 1; 
	//std::cout << "max_shift:" << max_bs << " bp_max:" << bp_max << std::endl; 
	size_t i; 
	uint32_t l1 = buf_.getlen();
	// bp: current bit plane
	for (bp = max_bs; bp > 0; bp--) {
		for (i = 0; i < genoUncVecSize_; i++) {
			GenoUnc &gu = genoUncVec_[i];
			if (gu.bp == 0) 
				continue;
			if ( gu.bh) {
				if (max_bs == bp) nRemaining--;
				continue;  
			}

			if ((gu.bs + gu.bp) == bp) {
				// gu.print();
				const SAMRecord &sr = samRecordDeque_[gu.read_idx];			
				d = sr.cidx;			

				if (((q[i]>>(gu.bp))<<(gu.bp)) + bit_mask[gu.bp-1]  <= qualitySegMax_c[d][gu.seg_id]) // otherwise bi will be zero with probability 1
				{
					bi = (q[i] & bit_mask[gu.bp-1]) ? 1 : 0;

					int prev_bit = 0; 
					if (gu.qual_idx && (prev_bp <= gu.bp))
						prev_bit = (prev_q & bit_mask[gu.bp-1]) ? 2:1; 

					int ctx = select_ctx(q[i], d, gu.seg_id, gu.bp, prev_bit);
					buf_.encode_bit(ctx, bi);

					if (buf_.getlen()-l >= nMaxBytes) // stop encoding when reaching maximum bit rate
						break;
				}
				prev_bp = gu.bp;
				prev_q = q[i];  
				gu.bp--;  
				if (gu.bp == 0) 
				 	nRemaining--;
			}
		}

		// check bitstream size inside fine scale encoding loop
		if (buf_.getlen()-l >= nMaxBytes) // stop encoding when reaching maximum bit rate
			break;
		if (!nRemaining) break;
	}
	//printf("Number of quality values encoded %d remaining: %d\n", uncompressedQualSize_, nRemaining);

	buf_.flush();
	compressedQualSize_ = buf_.getlen()-l;
	if (nRemaining) compressedQualSize_ --; // we don't want a potentially incomplete byte in bit-stream at lossy mode 

	if (genoUncVec_) free(genoUncVec_);
	if (q) free(q); 
}

void SQCCodec::decodeMappedQual(uint8_t *encoded_bytes)
{
	uint16_t bp_max, d, prev_bp, prev_q = 0;
	uint32_t nRemaining, bi, nMaxBytes;
	int16_t bp, max_bit_shifting, eof = 0;
	size_t i, nSamRecords = samRecordDeque_.size();
	uint8_t qr;
	
	// printf("Total number of SAM records: %zu\n", nSamRecords);
	for (i = 0; i < nSamRecords; i++) {
		SAMRecord &sr = samRecordDeque_[i];
		sr.qual = std::string(sr.qual.length(), '\0');
	}
	nRemaining = uncompressedQualSize_;
	nMaxBytes = (opt_.sqc_bpq >= 0) ? (uint32_t)(floor(nRemaining * opt_.sqc_bpq / 8.0)) : 0xffffffff; // bytes per quality
	// std::cout << "nMaxBytes: " << nMaxBytes << std::endl;

	genoUncVec_ = (GenoUnc *)calloc(uncompressedQualSize_, sizeof(GenoUnc));
	if (!genoUncVec_)
		 hts_log_error("SQC: memory failure");
	
	int hd = readQualitySegMinMaxTable(encoded_bytes);
	encoded_bytes += hd;
	compressedQualSize_ -= hd; 

	CABAC_BitstreamFile cabac_stream; 
	cabac_stream.Init(encoded_bytes,compressedQualSize_); 
	CABAC_ArithmeticDecoder arithmeticDecoder(&cabac_stream);
  	arithmeticDecoder.start();

	bp_max = updateQualitySegBPMax(false);
	decodeQualBitshifting(&arithmeticDecoder, max_bit_shifting);

	max_bit_shifting += bp_max;
	prev_bp = max_bit_shifting+1; 
	//std::cout << "max_shift:" << max_bit_shifting << " bp_max:" << bp_max << std::endl; 
 
	for (bp = max_bit_shifting; bp > 0; bp--) {
		for (i = 0; i < genoUncVecSize_; i++) {
			GenoUnc &gu = genoUncVec_[i];
			if (gu.bp == 0) 
				continue;	
			if (gu.bh) {
				if (max_bit_shifting == bp)	nRemaining --;
				continue;
			}		
			if ((gu.bs + gu.bp) == bp) {
				SAMRecord &sr = samRecordDeque_[gu.read_idx];
				qr = sr.qual[gu.qual_idx];
				d = sr.cidx;
				if (((qr >> (gu.bp)) << (gu.bp)) + bit_mask[gu.bp-1]  <= qualitySegMax_c[d][gu.seg_id]) //otherwise bi will be zero with probability 1
				{
					int prev_bit = 0; 
					if (gu.qual_idx && (prev_bp <= gu.bp))
						prev_bit = (prev_q & bit_mask[gu.bp-1]) ? 2:1; 

					int ctx = select_ctx(qr, d, gu.seg_id,gu.bp, prev_bit);

					eof = arithmeticDecoder.decodeBin(bi, &ctxt[ctx]);
					if (eof < 0)
						break; 
					
					if (bi != 0) // decoded bit is 1
						qr |= bit_mask[gu.bp-1];

					sr.qual[gu.qual_idx] = qr;
				}
				prev_bp = gu.bp; 
				prev_q = qr; 
				gu.bp--;
				if (gu.bp == 0)
					nRemaining--;
			}
			if (nMaxBytes < cabac_stream.getLen()) /* early terminating decoding based on bpq option for debugging purpose */
				break; 
		}

		if (eof < 0)
			break;
		if (nMaxBytes < cabac_stream.getLen()) /* early terminating decoding based on bpq option for debugging purpose */
			break; 
	}

	//printf("Number of quality values decoded %d remaining: %d\n", uncompressedQualSize_, nRemaining);
	
	// incomplete decoding in lossy mode
	for (i = 0; i < genoUncVecSize_; i++) {
		GenoUnc &gu = genoUncVec_[i];
		//if (gu.bp == 0) continue;

		if (!gu.bh) {
			SAMRecord &sr = samRecordDeque_[gu.read_idx];	
			uint16_t ql,qh; 
			uint8_t d = sr.cidx;
			ql = sr.qual[gu.qual_idx];
			if (gu.bp < opt_.sqc_cbr_bp ) {  // we start using context based reconstruction offset from sqc_cbr_bp  	 		
				qh = ql + (1<<gu.bp) - 1; 			
				qh = (qh < qualitySegMax_c[d][gu.seg_id]) ? qh : qualitySegMax_c[d][gu.seg_id]; 
				ql = qualitySegValueMap_[d][gu.seg_id][ql];
				qh = qualitySegValueMap_[d][gu.seg_id][qh];
				sr.qual[gu.qual_idx] = (uint8_t)(3*qh+ql)/4; 
			}
			else { 
				int offset_ctx = (ql >> gu.bp) + SQCCodec::offset_table[gu.bp - 1]; 
				sr.qual[gu.qual_idx]  = (uint8_t)cbQualRecOffset_[d][gu.seg_id][offset_ctx]; 
				//std::cout << "d:" << d+0 << " seg:" << gu.seg_id << " b:" << offset << " val:" << sr.qual[gu.qual_idx]+0  << std::endl; 
			}
		}
	}

#if debugging_code
	std::cout << std::endl;
	for (i = 0; i < nSamRecords; i++) {
		std::string q(samRecordDeque_.at(i).qual);
		for (auto & d : q) {
			d += (SQCCodec::qualityValueOffset_);
			std::cout << (int)d << ",";
		}
		std::cout << std::endl;
		std::cout << q << std::endl;
	}
#endif
	if (genoUncVec_) free(genoUncVec_);
}

int16_t SQCCodec::encodeQualBitshifting(DataBuffer & buf)
{
	size_t nSamRecords = samRecordDeque_.size(), nSamPileups = samPileupDeque_.size();
	size_t i, j, cigarIdx, cigarLen, opLen, idx, pileupIdx;
	uint32_t k, gu_idx = 0;
	uint16_t seg_pos, seg_id;
	int16_t d;
	bool bInsert;

	// printf("Total number of SAM records: %zu\n", nSamRecords); // number of SAM records
	if (!nSamPileups) {
		for (i = 0; i < nSamRecords; i++) {
			const SAMRecord &sr = samRecordDeque_[i];
#if debugging_code
			if (sr.isMapped() == true)
				hts_log_error("SQC: mapped SAMRecords in unmapped blocks");
#endif
			seg_pos = 0;
			seg_id = 0;
			d = sr.cidx;
			for (idx = 0; idx < sr.qual.size(); idx++) {
				genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), 0, i, idx, seg_id);
				genoUncVec_[gu_idx].bs = 0; 
				gu_idx++;
				seg_pos++;
				if (seg_pos >= opt_.segsize) {
					seg_id++;
					seg_pos = 0;
				}
			}
			// continue;
		}
		genoUncVecSize_ = gu_idx;
		return 0;
	}

	std::vector< std::pair<uint32_t, uint32_t> > gIdxInsertBitShifting; // use to update bit shifting of quality values of inserted bases
	std::vector<uint32_t> gIdxInsert; // used to record bit shifting of quality values of inserted bases, for encoding into buf_
	std::vector<uint32_t> gIdxMismatch; // used to record bit shifting of quality values of mismatched bases, for encoding into buf_

	// ------------------------------------------------------------------------------------------
	for (i = 0; i < nSamRecords; i++) {
		const SAMRecord &sr = samRecordDeque_[i];

		d = sr.cidx;
		if (sr.isMapped() == false) {
			seg_pos = 0;
			seg_id = 0;
			
			for (idx = 0; idx < sr.qual.size(); idx++) {
				genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), 0, i, idx, seg_id);
				genoUncVec_[gu_idx].bs = 0; 
				gu_idx++;
				seg_pos++;
				if (seg_pos >= opt_.segsize) {
					seg_id++;
					seg_pos = 0;
				}
			}
			continue;
		}

		cigarIdx = 0;
	    cigarLen = sr.cigar.length();
	    opLen = 0;  // length of current CIGAR operation
	    pileupIdx = sr.posMin - samPileupDeque_.posMin();
	    bInsert = false;
	    idx = 0;
		
	    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
	        if (isdigit(sr.cigar[cigarIdx])) {
	            opLen = opLen*10 + (size_t)sr.cigar[cigarIdx] - (size_t)'0';
	            continue;
	        }

#if debugging_code
	        std::cout << opLen << sr.cigar[cigarIdx] << " : " << idx << ", " << pileupIdx << std::endl;
#endif
	        switch (sr.cigar[cigarIdx]) {
	        case 'M':
	        case '=':
	        case 'X':
	        	for (j = 0; j < opLen; j++) {
	        		// idx is the actual loci of current base/qual in current read (grow with insertion and soft clipping)
	        		// k is its index on pileup (does not grow with insertion and soft clipping)
					seg_id = (idx/opt_.segsize);
		        	SAMPileup *spu = samPileupDeque_[pileupIdx];
#if debugging_code
	        		if (seg_id < 0 || seg_id > maxSegID_)
						hts_log_error("SQC: seg_id overflow");
	        		if (!spu)
						hts_log_error("SQC: empty SAMPileup");

	            	// to remove when debugged: assert equal by two ways of access: either through read, or through samPileup
	            	if (spu->seq[sr.idxPileup[idx]] != 'I' && spu->seq[sr.idxPileup[idx]] != 'D') {
						if ((sr.seq[idx] != spu->seq[sr.idxPileup[idx]]) || (sr.qual[idx] != spu->qual[sr.idxPileup[idx]])) {
							printf("(%c, %d), (%c, %d)\n", sr.seq[idx], (int)(sr.qual[idx]), spu->seq[sr.idxPileup[idx]], (int)(spu->qual[sr.idxPileup[idx]]));
							hts_log_error("SQC: inconsistent SAMRecord access");
						}
					}					
#endif
					genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), 0, i, idx, seg_id);
					int depth = spu->seq.length(); 
					int s;  
					for(s = 0; s < depth; s++) {
						if (spu->seq[s] != '=') break; 
					}
					if (s < depth )
						genoUncVec_[gu_idx].bs = 6; 
					else
						genoUncVec_[gu_idx].bs = 0; 
	            	gu_idx++;
					idx++; pileupIdx++;
	            }

        		bInsert = false;
	            break;
	        case 'I':
	        	bInsert = true;

	        	k = gu_idx + opLen;
	        	for (j = 0; j < opLen; j++) {
					seg_id = (idx/opt_.segsize);
#if debugging_code
	        		if (seg_id < 0 || seg_id > maxSegID_)
						hts_log_error("SQC: seg_id overflow");
#endif
					genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), pileupIdx, i, idx, seg_id);
					genoUncVec_[gu_idx].bs = 6; 

	        		gIdxInsertBitShifting.push_back( std::make_pair(gu_idx, k) );
	        		gu_idx++;
	        		idx++;
#if debugging_code
	        		std::cout << "inserted : ";
	        		genoUncVec_[gu_idx].print();
#endif
	        	}
	        	break;
	        case 'S':
	        	for (j = 0; j < opLen; j++) {
					seg_id = (idx/opt_.segsize);
#if debugging_code
	        		if (seg_id < 0 || seg_id > maxSegID_)
						hts_log_error("SQC: seg_id overflow");
#endif
					genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), pileupIdx, i, idx, seg_id);
					genoUncVec_[gu_idx].bs = 0; 
	        		gu_idx++;
	        		idx++;
	        	}
	        	if (bInsert)
	        		bInsert = false;
 	            break;
	        case 'D':
	        case 'N':
	            if (bInsert)
	            	bInsert = false;
	            pileupIdx += opLen;
	            break;
	        case 'H':
	        case 'P':
	        	if (bInsert)
	        		bInsert = false;
	            break;  // these have been clipped
	        default:
				hts_log_error("SQC: bad CIGAR string");
	        }
	        opLen = 0;
	    }
	}

	genoUncVecSize_ = gu_idx;

	return 6; // hardcode max_bs to 6 bits
}

int16_t SQCCodec::decodeQualBitshifting(CABAC_ArithmeticDecoder* p, int16_t & max_bs)
{	
	uint16_t d;
	uint32_t idx, k, gu_idx = 0;
	size_t i, j, cigarIdx, cigarLen, opLen, pileupIdx, nSamRecords = samRecordDeque_.size();
	bool bInsert;

	std::vector< std::pair<uint32_t, uint32_t> > gIdxInsertBitShifting; // use to update bit shifting of quality values of inserted bases
	std::vector<uint32_t> gIdxInsert; // used to record bit shifting of quality values of inserted bases, for encoding into buf_
	std::vector<uint32_t> gIdxMismatch; // used to record bit shifting of quality values of mismatched bases, for encoding into buf_

	// ------------------------------------------------------------------------------------------
	// init genotype and uncertainty at all loci, iterate over loci
	uint16_t seg_id,seg_pos;
	for (i = 0; i < nSamRecords; i++) {
		const SAMRecord &sr = samRecordDeque_[i];

		if (sr.isMapped() == false) {
			seg_pos = 0;
			seg_id = 0;
			d = sr.cidx;
			for (idx = 0; idx < sr.qual.size(); idx++) {
#if debugging_code
				if (seg_id < 0 || seg_id > maxSegID_)
					hts_log_error("SQC: seg_id overflow");
#endif
				genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), 0, i, idx, seg_id);
				genoUncVec_[gu_idx].bs = 0; 
				gu_idx++;
				seg_pos++;
				if (seg_pos >= opt_.segsize) {
					seg_id++;
					seg_pos = 0;
				}
			}
			continue;
		}

	    cigarIdx = 0;
	    cigarLen = sr.cigar.length();
	    opLen = 0;  // length of current CIGAR operation
	    pileupIdx = sr.posMin - samPileupDeque_.posMin();
		bInsert = false;
	    idx = 0;
		d = sr.cidx;
		
	    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
	        if (isdigit(sr.cigar[cigarIdx])) {
	            opLen = opLen*10 + (size_t)sr.cigar[cigarIdx] - (size_t)'0';
	            continue;
	        }

	        switch (sr.cigar[cigarIdx]) {
	        case 'M':
	        case '=':
			case 'X':
				for (j = 0; j < opLen; j++) {
					SAMPileup *spu = samPileupDeque_[pileupIdx];
					seg_id = (idx/opt_.segsize);
#if debugging_code
					if (seg_id < 0 || seg_id > maxSegID_)
						hts_log_error("SQC: seg_id overflow");
#endif
					
					genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), pileupIdx, i, idx, seg_id);
					int depth = spu->seq.length(); 
					int s;  
					for(s = 0; s < depth; s++) {
						if (spu->seq[s] != '=') break; 
					}
					if (s < depth )
						genoUncVec_[gu_idx].bs = 6; 
					else 
						genoUncVec_[gu_idx].bs = 0; 
					gu_idx++;
	            	idx++; pileupIdx++;
	            }
				bInsert = false;
	            break;
	        case 'I':
				bInsert = true;
				k = gu_idx + opLen;
				for (j = 0; j < opLen; j++) {
					seg_id = (idx/opt_.segsize);
#if debugging_code
					if (seg_id < 0 || seg_id > maxSegID_)
						hts_log_error("SQC: seg_id overflow");
#endif
					genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), pileupIdx, i, idx, seg_id);
					genoUncVec_[gu_idx].bs = 6; 
					gIdxInsertBitShifting.push_back( std::make_pair(gu_idx, k) );
					gu_idx++;
					idx++;
				}
				break;
	        case 'S':
	        	for (j = 0; j < opLen; j++) {
					seg_id = (idx/opt_.segsize);
#if debugging_code
					if (seg_id < 0 || seg_id > maxSegID_)
						hts_log_error("SQC: seg_id overflow");
#endif
					genoUncVec_[gu_idx].set_value(0, qualitySegBPMax_[d][seg_id], 0, (sr.seq[idx] == 'N'), pileupIdx, i, idx, seg_id);
					genoUncVec_[gu_idx].bs = 0; 
	        		idx++;
	        		gu_idx++;
	        	}
	            if (bInsert)
	            	bInsert = false;
	            break;
	        case 'D':
	        case 'N':
	            pileupIdx += opLen;
				if (bInsert)
					bInsert = false;
	            break;
	        case 'H':
	        case 'P':
				if (bInsert)
					bInsert = false;
	            break;  // these have been clipped
	        default:
				hts_log_error("SQC: bad CIGAR string");
	        }
	        opLen = 0;
	    }
	}
	genoUncVecSize_ = gu_idx;

	max_bs = 6; 
	return 0; 
}

int32_t SQCCodec::select_ctx(uint8_t q, uint8_t dir, uint16_t segment, int16_t bp, uint8_t prev_bit)
{
	int32_t context; 
	
	context = context_offset_Q; 
	context += dir; 
	context += (bp-1)*4;
	context += prev_bit*24; //prev_bit*6*4; 
	context += (q >> bp)*72; //(q >> bp)*3*6*4; 
	context += segment*4608; //segment*64*3*6*4; 
	return context; 
}
}  // namespace sqc
