//
//  gfastats-input.h
//  gfastats
//
//  Created by Giulio Formenti on 1/16/22.
//

#ifndef GFASTATS_INPUT_H
#define GFASTATS_INPUT_H

class Input {
    
    UserInput userInput;
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    bool stopStream = false;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    StreamObj streamObj;
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    std::shared_ptr<std::istream> stream;
    
    std::vector<Instruction> instructions;
    
    unsigned int begin = 0, end = 0;
    
public:
    
    void load(UserInput userInput);
    
    void read(InSequences& inSequence);
    
    void read(InReads& inReads);
    
};

#endif /* GFASTATS_INPUT_H */
