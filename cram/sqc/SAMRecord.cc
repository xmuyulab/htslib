/** @file SAMRecord.cc
 *  @brief This file contains the implementation of the SAMRecord class.
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

// This code is modified from CALQ developed by Jan Voges
// https://github.com/voges/calq


#include <string.h>

#include "SAMRecord.h"
#include "htslib/hts_log.h"

namespace sqc {

SAMRecord::SAMRecord(char *fields[NUM_FIELDS])
    : flag((uint16_t)atoi(fields[1])),
      pos((uint32_t)atoi(fields[3])),
      mapq((uint8_t)atoi(fields[4])),
      cigar(fields[5]),
      seq(fields[9]),
      qual(fields[10]),
      posMin(0),
      posMax(0),
	  qualMin(128),
	  qualMax(0),
      mapped_(false) {
    check();

    if (mapped_ == true) {
        // Compute 0-based first position and 0-based last position this record
        // is mapped to on the reference used for alignment
        posMin = pos - 1;
        posMax = pos - 1;

        size_t cigarIdx = 0;
        size_t cigarLen = cigar.length();
        uint32_t opLen = 0;  // length of current CIGAR operation

        for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
            if (isdigit(cigar[cigarIdx])) {
                opLen = opLen * 10 + (uint32_t)cigar[cigarIdx] - (uint32_t)'0';
                continue;
            }
            switch (cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                posMax += opLen;
                break;
            case 'I':
            case 'S':
                break;
            case 'D':
            case 'N':
                posMax += opLen;
                break;
            case 'H':
            case 'P':
                break;  // these have been clipped
            default:
                hts_log_error("SQC: bad CIGAR string");
            }
            opLen = 0;
        }

        posMax -= 1;
    }

    for (auto c : qual) {
    	if (c > qualMax) qualMax = c;
    	if (c < qualMin) qualMin = c;
    }

    if (mapped_ && (posMax <= posMin)) {
    	printf("Read length not positive, (posMin, posMax) = (%u, %u), cigar = %s\n", posMin, posMax, std::string(cigar).c_str());        
    }
    readLen = posMax - posMin;
    ref_id = -1;
}

SAMRecord::SAMRecord(int32_t flags_, int32_t ref_id_, int32_t apos_, int32_t aend_, int32_t len_, const char * seq_, const char * ext_cig_)
	: flag(flags_),
	  readLen(len_),
	  ref_id(ref_id_),
	  posMin(apos_),
      posMax(aend_),
	  qualMin(128),
	  qualMax(0),
      mapped_(false) {

	seq = std::string(seq_, readLen);
	cigar = std::string(ext_cig_);
	qual = std::string(seq_, readLen);
    pos = apos_;

	check();
}

SAMRecord::~SAMRecord(void) {}

void SAMRecord::addToPileupQueue(SAMPileupDeque &samPileupDeque_, const std::string &qualPred) {
    if (samPileupDeque_.empty() == true) 
        hts_log_error("SQC: samPileupDeque is empty");
    if ((samPileupDeque_.posMin() > posMin) || (samPileupDeque_.posMax() < posMax))
        hts_log_error("SQC: record not in the range of samPileupDeque");

    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0;  // length of current CIGAR operation
    size_t i, idx = 0; // actual index of seq and qual in current read, grow with insert and softclipping
    size_t pileupIdx = posMin - samPileupDeque_.posMin();
    bool bInsert = false; // true if previous cigar is I(nsert)

    idxPileup.clear();
    idxPileup.resize(qual.size(), 0);

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen*10 + (size_t)cigar[cigarIdx] - (size_t)'0';
            continue;
        }

        switch (cigar[cigarIdx]) {
        case 'M':
        case '=':
        case 'X':
        	// for the first base
        	if (!samPileupDeque_[pileupIdx]) {
        		samPileupDeque_.pileups_[pileupIdx] = new SAMPileup;
#if debugging_code
        		if (!samPileupDeque_[pileupIdx])
                    hts_log_error("SQC: memory failure");
#endif
        	}
        	if (bInsert) // previous cigar is Insert, combine with current base
        		samPileupDeque_[pileupIdx]->setValues(samPileupDeque_.posMin() + pileupIdx, 'I', qualPred[idx]);
        	else
                samPileupDeque_[pileupIdx]->setValues(samPileupDeque_.posMin() + pileupIdx, cigar[cigarIdx], qualPred[idx]);

        	idxPileup[idx] = samPileupDeque_[pileupIdx]->size() - 1;
        	idx++; pileupIdx++;

            for (i = 1; i < opLen; i++) {
            	if (!samPileupDeque_[pileupIdx]) {
            		samPileupDeque_.pileups_[pileupIdx] = new SAMPileup;
            		if (!samPileupDeque_[pileupIdx])
                        hts_log_error("SQC: memory failure");
            	}            	
                samPileupDeque_[pileupIdx]->setValues(samPileupDeque_.posMin() + pileupIdx, cigar[cigarIdx], qualPred[idx]);
            	idxPileup[idx] = samPileupDeque_[pileupIdx]->size() - 1;
                idx++; pileupIdx++;
            }
            bInsert = false;
            break;
        case 'I':
        	bInsert = true;
        	idx += opLen;
        	break;
        case 'S': // soft clipping quals are not added to pileups
        	bInsert = false;
            idx += opLen;
            break;
        case 'D':
        	if (!samPileupDeque_[pileupIdx]) {
        		samPileupDeque_.pileups_[pileupIdx] = new SAMPileup;
        		if (!samPileupDeque_[pileupIdx])
                    hts_log_error("SQC: memory failure");
        	}
        	samPileupDeque_[pileupIdx]->setValues(samPileupDeque_.posMin() + pileupIdx, 'D', qualPred[idx]);
        case 'N':
            pileupIdx += opLen;
            bInsert = false;
            break;
        case 'H':
        case 'P':
        	bInsert = false;
            break;  // these have been clipped
        default:
            hts_log_error("SQC: bad CIGAR string");
        }

        opLen = 0;
    }
}

bool SAMRecord::isMapped(void) const {
    return mapped_;
}

void SAMRecord::printLong(void) const {
    printShort();
    printf("isMapped: %d, ", mapped_);
    printf("posMin: %d, ", posMin);
    printf("posMax: %d, ", posMax);
    printf("readLen: %d, ", readLen);
    printf("(seq, qual) pairs in samPileup: %zu\n", idxPileup.size());
    for (size_t i = 0; i < idxPileup.size(); i++)
    	printf("%zu,", idxPileup[i]);
    printf("\n");
    printf("%s\n", qual.c_str());
}

void SAMRecord::printShort(void) const {
    printf("%d\t", flag);
    printf("%d\t", pos);
    printf("%d\t", mapq);
    printf("%s\t", cigar.c_str());
    printf("%s\t", seq.c_str());
    printf("%s\t", qual.c_str());
    printf("%s\t", opt.c_str());
    printf("\n");
}

void SAMRecord::check(void) {
    if (cigar.empty() == true)
        hts_log_error("SQC: empty CIGAR string");
    if (seq.empty() == true) 
        hts_log_error("SQC: empty SEQ string");
    if (qual.empty() == true) 
        hts_log_error("SQC: empty QUAL string");

    // Check if this record is mapped
    if ((flag & 0x4) != 0) {
        mapped_ = false;
    } else {
        mapped_ = true;
        if (cigar == "*" || seq == "*" || qual == "*") 
            hts_log_error("SQC: corrupted SAMRecord");
    }
    cidx =  ((((flag & 0x10)) ? 2 : 0) + (((flag & 0x40)) != 0));
}
}  // namespace sqc

