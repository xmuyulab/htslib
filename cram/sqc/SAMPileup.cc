/** @file SAMPileup.cc
 *  @brief This file contains the implementation of the SAMPileup class.
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

// This code is modified from CALQ developed by Jan Voges
// https://github.com/voges/calq


#include "SAMPileup.h"

#include <stdio.h>

namespace sqc {

SAMPileup::SAMPileup(void) : pos(0), seq(""), qual("")
{
	score.clear(); 
}

SAMPileup::~SAMPileup(void) {}

bool SAMPileup::empty(void) const {
    return seq.empty();
}

void SAMPileup::clear(void) {
    pos = 0;
    seq = "";
    qual = "";
	score.clear(); 
}

void SAMPileup::print(void) const {
	printf("(pos: %6d, len: %3zu) %s %s\n", pos, seq.size(), seq.c_str(), qual.c_str());
}

void SAMPileup::printSeq(void) const {
    printf("%6d: %s\n", pos, seq.c_str());
}

void SAMPileup::setValues(uint32_t p, char b, char q) {
	pos = p;
	seq += b;
	qual += q;
}
}  // namespace sqc
