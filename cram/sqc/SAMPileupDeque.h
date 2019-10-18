/** @file SAMPileupDeque.h
 *  @brief This file contains the definition of the SAMPileupDeque class.
 *
 *  Copyright (c) 2018-2019, Aginome Scientific
 *  All rights reserved.
 */

// This code is modified from CALQ developed by Jan Voges
// https://github.com/voges/calq


#ifndef SQC_SAMPILEUPDEQUE_H_
#define SQC_SAMPILEUPDEQUE_H_

#include <deque>

#include "SAMPileup.h"

namespace sqc {

class SAMPileupDeque {
    friend class SAMRecord;

 public:
    SAMPileupDeque(void);
    ~SAMPileupDeque(void);

    void clear(void);
    bool empty(void) const;
    size_t length(void) const;
    SAMPileup * operator[](const size_t &n) const;
    size_t size(void) const;
    void print(void) const;

    uint32_t posMax(void) const;
    uint32_t posMin(void) const;

    void setPosMax(uint32_t posMax);
    void setPosMin(uint32_t posMin);

 private:
    std::vector<SAMPileup *> pileups_;
    uint32_t posMax_;
    uint32_t posMin_;
};

}  // namespace sqc

#endif  // SQC_SAMPILEUPDEQUE_H_

