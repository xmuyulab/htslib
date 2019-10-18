/** @file sqc_wrapper.h
 *  @brief This file contains the definition of the sqc_wrapper functions
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

#ifndef _SQC_WRAPPER_H_
#define _SQC_WRAPPER_H_

#include "cram/cram_samtools.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cram_block cram_block;
typedef struct cram_fd cram_fd;
typedef struct cram_record cram_record;

// SQC options for SQC scalable quality value compression
typedef struct sqc_option {
	int sqc_enc;
	int sqc_dec;
	int sqc_cut; // bit stream truncation
	double sqc_bpq;
	unsigned int sqc_cbr_bp; // lowest bitplane that we use context based decoding offset 
	unsigned int segsize; // length of a segment in a read/quality string
} sqc_option;

// wrapper for encoder
struct sqc_SQCCodec {
	void *obj;
};
typedef struct sqc_SQCCodec sqc_SQCCodec_t;

// constructor and deconstructor
sqc_SQCCodec_t* sqc_SQCCodec_create();
void sqc_SQCCodec_destroy(sqc_SQCCodec_t *e);

// public member functions
void sqc_updateQualitySegMinMax(sqc_SQCCodec_t *e, const bam_seq_t *b, sqc_option opt);
void sqc_addRecordToBlock(sqc_SQCCodec_t *e, const bam_hdr_t *hdr, const bam_seq_t *b, const char * ext_cig, sqc_option opt);
void sqc_addRecordToBlock_cr(sqc_SQCCodec_t *e, const cram_record * cr, const char * seq, const char * ext_cig, sqc_option opt);

int sqc_finishBlock(sqc_SQCCodec_t *e, cram_block *b, sqc_option opt);
int sqc_decodeBlock(sqc_SQCCodec_t *e, cram_block *b, cram_block *bout, sqc_option opt);
int sqc_cutBlock(sqc_SQCCodec_t *e, cram_block *b, double bpq_out);
int sqc_getTruncatedBitstream(const sqc_SQCCodec_t *e, uint8_t *b);
int sqc_setTruncatedBitstream(sqc_SQCCodec_t *e, const uint8_t *b, int len);

void sqc_initRefOffsetVec(sqc_SQCCodec_t *e, int nref);

#ifdef __cplusplus
}
#endif

#endif
