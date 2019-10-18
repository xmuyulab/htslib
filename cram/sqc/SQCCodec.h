/** @file SQCCodec.h
 *  @brief This file contains the definition of the SQCCodec class.
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

#ifndef SQC_CODEC_H_
#define SQC_CODEC_H_

#include <chrono>
#include <deque>
#include <map>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <math.h>

#include "SAMPileupDeque.h"
#include "SAMRecord.h"
#include "CABAC_ArithmeticEncoder.h"
#include "CABAC_ArithmeticDecoder.h"
#include "CABAC_BitstreamFile.h"
#include "ContextModel.h"

#include "cram/cram_structs.h"

namespace sqc {

#define DATA_BUF_SIZE 1000000 // 1MB
#define MAX_SEGMENTS 100 // maximum number of segments per read/quality string
#define MAX_NQVALUES 94 // maximum number of valid quality values, Sanger format [0,93]

typedef std::vector < std::vector<uint32_t> > IndexList;

struct GenoUnc {
	double score;

	uint8_t bp; // bp: bit plane
    uint8_t bs;	// bs: bit shifting
	uint8_t bh; // bh: whether corresponding base is 'N'
	uint32_t pile_idx;
	uint32_t read_idx;
	uint32_t qual_idx;
	uint8_t seg_id; // segment id

	GenoUnc() : score(0), bp(0), bs(0), bh(0), pile_idx(0), read_idx(0), qual_idx(0), seg_id(0)
	{}

    GenoUnc(double score_, int16_t bp_, int16_t bs_, int16_t bh_, uint32_t pile_idx_, uint32_t read_idx_, uint32_t qual_idx_, uint16_t seg_id_) :
    score (score_), bp(bp_), bs(bs_), bh(bh_), pile_idx(pile_idx_), read_idx(read_idx_), qual_idx(qual_idx_), seg_id(seg_id_)
    {}

    GenoUnc(const GenoUnc & g) // copy constructor
    { score = g.score; bp = g.bp; bs = g.bs; bh = g.bh; pile_idx = g.pile_idx; read_idx = g.read_idx; qual_idx = g.qual_idx; seg_id = g.seg_id; }

    void set_value(double score_, uint8_t bp_, uint8_t bs_, uint8_t bh_, uint32_t pile_idx_, uint32_t read_idx_, uint32_t qual_idx_, uint8_t seg_id_)
    {
        score = score_;
        bp = bp_;
        bs = bs_;
        bh = bh_;
        pile_idx = pile_idx_;
        read_idx = read_idx_;
        qual_idx = qual_idx_;
        seg_id = seg_id_;
    }
	void print() { std::cout << pile_idx << "," << read_idx << "," << qual_idx << "," << bp << "," << bs << "," << bh << "," << score << std::endl; }
};

class DataBuffer {
public:
	DataBuffer() { dat = 0; ptr = 0; size = 0; total_byte = 0; }
	~DataBuffer() { if (dat) free(dat); }

	uint8_t *dat; // encoded bitstream containing all mapped, separately coded, and unmapped quality values
    uint8_t *ptr; // pointer to current position in data

private:
	size_t size; // allocated size of buffer data
    uint32_t total_byte; 

    ContextModel *ctxt; 
    CABAC_BitstreamFile cabac_stream; 
    CABAC_ArithmeticEncoder arithmeticEncoder; 

public:
	int init() {
        dat = (uint8_t *)calloc(DATA_BUF_SIZE, sizeof(uint8_t));
		if (!dat) return -1;

		ptr = dat;
		size = DATA_BUF_SIZE;
        total_byte = 0;
    	return 0;
	}
    int init_ac(ContextModel *context) {
        ctxt = context; 
        cabac_stream.Init(ptr,0); 
        arithmeticEncoder.setBitstream(&cabac_stream); 
        arithmeticEncoder.start();
        return 0;
    }

    int encode_bit(int cx, int bit)
    {
        arithmeticEncoder.encodeBin(bit, &ctxt[cx]); 

        total_byte = ptr-dat+cabac_stream.getLen(); 
        if (total_byte  >= size){
            uint8_t * dat0;
			dat0 = (uint8_t *)realloc(dat, size + DATA_BUF_SIZE);  
			if (!dat0) return -1;
			ptr = dat0 + (ptr - dat) ; 
            dat = dat0; 
            size += DATA_BUF_SIZE;
            cabac_stream.reAlloc(ptr); 
		}
        return 0; 
    }    
    unsigned int flush()
    {
        int ac_len; 
        arithmeticEncoder.finish();
        ac_len = cabac_stream.getLen(); 
        total_byte = ptr-dat+ac_len; 
        ptr += ac_len; 
        return ac_len; 
    }
    uint32_t getlen() const
    {
        return total_byte; 
    }
    // direct read and write for bitstream truncation
    int setBuf(const uint8_t *p, uint32_t l) 
    {
        if (size == 0) // buffer not initialized 
            init();
        if (l >= size) {
            uint8_t * dat0;
			dat0 = (uint8_t *)realloc(dat, size + DATA_BUF_SIZE);  
			if (!dat0) return -1;
			ptr = dat0 + (ptr - dat) ; 
            dat = dat0; 
            size += DATA_BUF_SIZE;
		}
        memcpy(ptr, p, l * sizeof(uint8_t));
        total_byte = l;
        // for (int i = 0; i < l; i++)
        //     printf("%d, ", (uint16_t)p[i]);
        // printf("\n");
        return total_byte;
    }
    int getBuf(uint8_t *p, uint32_t l) const
    {
        if (l >= total_byte) {
            std::cout << "warning: buffer size smaller than requested" << std::endl;
            l = total_byte;
        }
        memcpy(p, ptr, l * sizeof(uint8_t));
        // for (int i = 0; i < l; i++)
        //     printf("%d, ", (uint16_t)ptr[i]);
        // printf("\n");
        return l;
    }
    int getBuf(uint8_t *p) const
    {
        memcpy(p, ptr, total_byte * sizeof(uint8_t));
        // for (int i = 0; i < total_byte; i++)
        //     printf("%d, ", (uint16_t)ptr[i]);
        // printf("\n");
        return total_byte;
    }
};

class SQCCodec {
public:
    explicit SQCCodec();
    ~SQCCodec(void);

    void updateQualitySegMinMax(const bam_seq_t *b, sqc_option opt);
    void addUnmappedRecordToBlock(const SAMRecord &samRecord, sqc_option opt);
    void addMappedRecordToBlock(SAMRecord &samRecord, sqc_option opt);

    void encodeBlock(cram_block *b);
    void decodeBlock(cram_block *b, cram_block *bout);
    size_t writeBlock(cram_block *b);

    void truncateBitstream(cram_block *b, double out_bpq);
    int getTruncatedBitstream(uint8_t *b) const; // from sqc buffer to b
    int setTruncatedBitstream(const uint8_t *b, int len); // from b to sqc buffer

    uint32_t nrMappedRecords(void) const { return nrMappedRecords_; }
    uint32_t nrUnmappedRecords(void) const { return nrUnmappedRecords_; }
    uint32_t nrRecords(void) const { return (nrMappedRecords_ + nrUnmappedRecords_); }

    void initRefOffsetVec(int nref);

    sqc_option opt_;
private:
    // Data buffer for compressed byte string
    DataBuffer buf_;
    
    uint32_t nrMappedRecords_;
    uint32_t nrUnmappedRecords_;
    uint32_t uncompressedQualSize_;
    uint32_t compressedQualSize_;
    uint32_t headerSize_;

    unsigned int qualitySegMin_[4][MAX_SEGMENTS] = { }; // per segment value
    unsigned int qualitySegMax_[4][MAX_SEGMENTS] = { };
    unsigned int qualitySegMax_c[4][MAX_SEGMENTS] = { }; // max quality after range compression
    double cbQualRecOffset_[4][MAX_SEGMENTS][64] = {}; // context based decoding reconstruction offset
    unsigned int maxSegID_;
    unsigned int qualitySegBPMax_[4][MAX_SEGMENTS] = { };
    unsigned int qualitySegValueValid_[4][MAX_SEGMENTS][MAX_NQVALUES] = { };
    unsigned int qualitySegValueMap_[4][MAX_SEGMENTS][MAX_NQVALUES] = { };
    const int offset_table[6] = {0,32,32+16,32+16+8,32+16+8+4, 32+16+8+4+2} ;
    
    static int qualityValueOffset_;

    SAMPileupDeque samPileupDeque_;
    std::vector<SAMRecord> samRecordDeque_;

    GenoUnc * genoUncVec_;
    size_t genoUncVecSize_;

    const uint32_t context_offset_I = 0; 
    const uint32_t context_offset_X = 8; 
    const uint32_t context_offset_Q = 16; 
    ContextModel ctxt[64*3*6*4*10+16];

private:
    void encodeMappedQual();
    void decodeMappedQual(uint8_t *p);

    int16_t encodeQualBitshifting(DataBuffer & buf);
    int16_t writeQualitySegMinMaxTable(uint8_t * buf);
    uint16_t updateQualitySegBPMax(bool enc);
    int16_t decodeQualBitshifting(CABAC_ArithmeticDecoder* p, int16_t & max_bs); 
    int16_t readQualitySegMinMaxTable(uint8_t *p);
    int32_t select_ctx(uint8_t q, uint8_t dir, uint16_t segment, int16_t bp, uint8_t prev_bit);

    uint16_t maxBitplane(uint16_t minval, uint16_t maxval)
    {
        int bp = 0; 
        while (minval+(1<<bp)-1 < maxval)
            bp++; 
        return bp;  
    };

    std::vector<int32_t> refOffsetVec_; // for multi-seq
};
}  // namespace sqc

#endif  // SQC_CODEC_H_

