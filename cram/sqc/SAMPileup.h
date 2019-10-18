/** @file SAMPileup.h
 *  @brief This file contains the definition of the SAMPileup class.
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

// This code is modified from CALQ developed by Jan Voges
// https://github.com/voges/calq

#ifndef SQC_SAMPILEUP_H_
#define SQC_SAMPILEUP_H_

#include <string>
#include <vector>

namespace sqc {

class SAMPileup {
 public:
    SAMPileup(void);
    ~SAMPileup(void);

    bool empty(void) const;
    void clear(void);
    void print(void) const;
    void printSeq(void) const;

    void setValues(uint32_t pos, char base, char qual);

    size_t size(void) const {
    	return seq.size();
    }

    uint32_t pos;  // 0-based position of this pileup
    std::string seq;
    std::string qual;
    std::vector<double> score; 
};

}  // namespace sqc

#endif  // SQC_SAMPILEUP_H_

