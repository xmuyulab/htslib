/** @file sqc_wrapper.cc
 *  @brief This file contains the implementation of the sqc_wrapper functions
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

#include "sqc_wrapper.h"
#include "SAMRecord.h"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "SQCCodec.h"

#include "htslib/thread_pool.h"

using namespace sqc;

static void parseLine(char *fields[SAMRecord::NUM_FIELDS], char *line) {
    char *c = line;
    char *pc = c;
    int f = 0;

    // for mandatory fields (both mapped and unmapped records)
    while (f <= 10) {
    	if (*c == '\t' || *c == '\0') {
    		*c = '\0';
    		fields[f] = pc;
    		f++;
    		pc = c + 1;
    	}
    	c++;
    }
}

sqc_SQCCodec_t* sqc_SQCCodec_create()
{
	sqc_SQCCodec_t *e;
	SQCCodec *obj;

	e = (__typeof__(e))malloc(sizeof(sqc_SQCCodec_t));
	if (!e) return NULL;
	obj = new SQCCodec();
	if (!obj) return NULL;
	e->obj = obj;

	return e;
}

void sqc_SQCCodec_destroy(sqc_SQCCodec_t *e)
{
	if (e == NULL) return;
	delete static_cast<SQCCodec *>(e->obj);
	free(e);
}

void sqc_updateQualitySegMinMax(sqc_SQCCodec_t *e, const bam_seq_t *b, sqc_option opt)
{
	if (e == NULL) return;

	SQCCodec *obj;
	obj = static_cast<SQCCodec *>(e->obj);
	obj->updateQualitySegMinMax(b, opt);
}

void sqc_initRefOffsetVec(sqc_SQCCodec_t *e, int nref)
{
	SQCCodec *obj;
	obj = static_cast<SQCCodec *>(e->obj);
	obj->initRefOffsetVec(nref);
}

void sqc_addRecordToBlock_cr(sqc_SQCCodec_t *e, const cram_record *cr, const char *seq, const char *ext_cig, sqc_option opt)
{
	if (e == NULL) return;

	SQCCodec *obj = static_cast<SQCCodec *>(e->obj);
	SAMRecord samRecord(cr->flags, cr->ref_id, (cr->apos-1), (cr->aend-1), cr->len, seq, ext_cig); 
	 if ((cr->flags & 0x4) != 0) {
	 	obj->addUnmappedRecordToBlock(samRecord, opt);
	 	return;
	 }
	 obj->addMappedRecordToBlock(samRecord, opt);
}

void sqc_addRecordToBlock(sqc_SQCCodec_t *e, const bam_hdr_t *hdr, const bam_seq_t *b, const char *ext_cig, sqc_option opt)
{
	if (e == NULL) return;

	SQCCodec *obj;
	char *fields[SAMRecord::NUM_FIELDS];

	kstring_t ks = { 0, 0, NULL };
	if (sam_format1(hdr, b, &ks) < 0) {
		hts_log_error("SQC: error formating SAMRecord");
		return;
	}
	parseLine(fields, ks.s);

	SAMRecord samRecord(fields);
	samRecord.ref_id = b->core.tid;

	// std::cout << samRecord.posMin << ", " << samRecord.posMax << std::endl;
//	samRecord.printShort();
	obj = static_cast<SQCCodec *>(e->obj);

	if (((b->core.flag) & 0x4) != 0)
		obj->addUnmappedRecordToBlock(samRecord, opt);
	else {
		if (ext_cig) {
			if (!opt.sqc_enc || opt.sqc_dec)
				hts_log_error("SQC: extended CIGAR string non available");

			samRecord.cigar = std::string(ext_cig);
		} else {
			hts_log_error("SQC: extended CIGAR string non available");
		}

		// std::cout << samRecord.cigar << std::endl;

		obj->addMappedRecordToBlock(samRecord, opt);
	}
	if (ks.s) free(ks.s);
}

int sqc_finishBlock(sqc_SQCCodec_t *e, cram_block *b, sqc_option opt)
{
	if (e == NULL) return -1;
		
	if (opt.sqc_cut) {
		b->method = AGI_SQC;
		b->comp_size = sqc_getTruncatedBitstream(e, b->data);
		return 0;
	}

	SQCCodec *obj;
	obj = static_cast<SQCCodec *>(e->obj);
	obj->opt_ = opt;
	obj->encodeBlock(b);
	return 0;
}

//int sqc_decodeBlock(sqc_SQCCodec_t *e, cram_block *b, bam1_t **bbuf, sqc_option opt)
int sqc_decodeBlock(sqc_SQCCodec_t *e, cram_block *b, cram_block *bout, sqc_option opt)
{
	if (e == NULL) return -1;

	SQCCodec *obj;
	obj = static_cast<SQCCodec *>(e->obj);
	obj->opt_ = opt;
	obj->decodeBlock(b, bout);
	return 0;
}

int sqc_cutBlock(sqc_SQCCodec_t *e, cram_block *b, double bpq_out)
{
	if (e == NULL) return -1;
	SQCCodec *obj = static_cast<SQCCodec *>(e->obj);
	obj->truncateBitstream(b, bpq_out);
	return 0;
}

int sqc_getTruncatedBitstream(const sqc_SQCCodec_t *e, uint8_t *b)
{
	if (e == NULL) return -1;
	SQCCodec *obj = static_cast<SQCCodec *>(e->obj);
	return obj->getTruncatedBitstream(b);
}

int sqc_setTruncatedBitstream(sqc_SQCCodec_t *e, const uint8_t *b, int len)
{
	if (e == NULL) return -1;
	SQCCodec *obj = static_cast<SQCCodec *>(e->obj);
	return obj->setTruncatedBitstream(b, len);
}