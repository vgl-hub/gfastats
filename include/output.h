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
    bool seqReport(InSequences &inSequences, int &outSequence_flag);
    
    bool outFile(InSequences &inSequences, int splitLength, std::string &outSeq);
    
    bool outSize(InSequences &inSequences, char &sizeOutType);
    
    bool outCoord(InSequences &inSequences, char bedOutType);
    
    bool reportStats(InSequences &inSequences, unsigned long long int gSize, InReads& inReads);
    
    bool nstarReport(InSequences &inSequences, unsigned long long int gSize);
    
};

#endif /* OUTPUT_H */
