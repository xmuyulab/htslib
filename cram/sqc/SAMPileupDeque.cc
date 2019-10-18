/** @file SAMPileupDeque.cc
 *  @brief This file contains the implementation of the SAMPileupDeque class.
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

// This code is modified from CALQ developed by Jan Voges
// https://github.com/voges/calq


#include "SAMPileupDeque.h"
#include "htslib/hts_log.h"

namespace sqc {

SAMPileupDeque::SAMPileupDeque(void) : pileups_() , posMax_(0) , posMin_(0) { pileups_.clear();}

SAMPileupDeque::~SAMPileupDeque(void)
{
	if (!pileups_.empty()) {
		for (auto & g : pileups_) {
			if (g) delete g;
		}
	}
}

void SAMPileupDeque::clear(void) {
    pileups_.clear();
    posMax_ = 0;
    posMin_ = 0;
}

bool SAMPileupDeque::empty(void) const {
	return pileups_.empty();
}

size_t SAMPileupDeque::length(void) const {
    return posMax_ - posMin_ + 1;
}

SAMPileup * SAMPileupDeque::operator[](const size_t &n) const {
    return pileups_[n];
}

size_t SAMPileupDeque::size(void) const {
    return pileups_.size();
}

void SAMPileupDeque::print(void) const {
    if (pileups_.empty()) {
        hts_log_error("SQC: deque is empty");
    }

   	for (auto const &samPileup : pileups_) {
   		if (samPileup) samPileup->print();
   	}
}

uint32_t SAMPileupDeque::posMax(void) const {
    return posMax_;
}

uint32_t SAMPileupDeque::posMin(void) const {
    return posMin_;
}

void SAMPileupDeque::setPosMax(uint32_t posMax) {
    if (posMax < posMax_) {
        hts_log_error("SQC: posMax out of range");
    }
    posMax_ = posMax;
    pileups_.resize(length());
}

void SAMPileupDeque::setPosMin(uint32_t posMin) {
    if (posMin < posMin_) {
        hts_log_error("SQC: posMin out of range");
    }

    if (pileups_.empty()) {
        posMin_ = posMin;
    } else {
        hts_log_error("SQC: SAMPileupDeque is non empty");
    }
}
}  // namespace sqc

