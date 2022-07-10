//
//  gfastats-output.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFASTATS_OUTPUT_H
#define GFASTATS_OUTPUT_H

//classes
class Report {

private:
    unsigned int counter = 0;
    
public:
    bool seqReport(InSequences &inSequences, InSegment &inSegment, int &outSequence_flag);
    
    bool outFile(InSequences &inSequences, InSegment &inSegment, int splitLength, std::string &outSeq);
    
    bool outSize(InSequences &inSequences, InSegment &inSegment, char &sizeOutType);
    
    bool outCoord(InSequences &inSequences, InSegment &inSegment, char bedOutType);
    
    bool reportStats(InSequences &inSequences, unsigned long long int gSize, int bedOutType);
    
    bool nstarReport(InSequences &inSequences, unsigned long long int gSize);
    
};

#endif /* GFASTATS_OUTPUT_H */
