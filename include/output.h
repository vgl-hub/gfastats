//
//  gfastats-output.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef OUTPUT_H
#define OUTPUT_H

//classes
class Report {

private:
    unsigned int counter = 0;
    
public:
    bool seqReport(InSequences &inSequences, InSegment &inSegment, int &outSequence_flag);
    
    bool outFile(InSequences &inSequences, InSegment &inSegment, int splitLength, std::string &outSeq);
    
    bool outSize(InSequences &inSequences, InSegment &inSegment, char &sizeOutType);
    
    bool outCoord(InSequences &inSequences, InSegment &inSegment, char bedOutType);
    
    bool reportStats(InSequences &inSequences, unsigned long long int gSize, int bedOutType, InReads& inReads);
    
    bool nstarReport(InSequences &inSequences, unsigned long long int gSize);
    
};

#endif /* OUTPUT_H */
