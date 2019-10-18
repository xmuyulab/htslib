/** @file SAMRecord.h
 *  @brief This file contains the definition of the SAMRecord class.
 */

// This code is modified from CALQ developed by Jan Voges
// https://github.com/voges/calq


#ifndef SQC_SAMRECORD_H_
#define SQC_SAMRECORD_H_

#include <inttypes.h>

#include <deque>
#include <string>
#include <vector>

#include "SAMPileupDeque.h"

namespace sqc {

class SAMRecord {
 public:
    static const int NUM_FIELDS = 12;

    SAMRecord(char *fields[NUM_FIELDS]);
    SAMRecord(int32_t flags_, int32_t ref_id_, int32_t apos_, int32_t aend_, int32_t len_, const char * seq_, const char * ext_cig_);
    ~SAMRecord(void);

    void addToPileupQueue(SAMPileupDeque &samPileupDeque, const std::string &qualPred);

    bool isMapped(void) const;
    void printLong(void) const;
    void printShort(void) const;

    uint16_t    flag;   // bitwise FLAG (uint16_t)   
    uint32_t    pos;    // 1-based leftmost mapping POSition (uint32_t)    
    uint8_t     mapq;   // MAPping Quality (uint8_t)
    std::string cigar;  // CIGAR string
    std::string seq;    // segment SEQuence
    std::string qual;   // QUALity scores
    std::string opt;    // OPTional information

    uint32_t posMin;  // 0-based leftmost mapping position
    uint32_t posMax;  // 0-based rightmost mapping position

    uint16_t qualMin;
    uint16_t qualMax;
    uint16_t readLen;
    uint16_t ref_id; // ref_id is bam record core->tid

    uint8_t cidx; // context index calculat from direction and first/second read

    std::vector<size_t> idxPileup;

 private:
    void check(void);
    SAMPileup * setSAMPileup(SAMPileup *spu, uint32_t pos, char base, char qual, uint16_t rpos);

 private:
    bool mapped_;
};

}  // namespace sqc

#endif  // SQC_SAMRECORD_H_

