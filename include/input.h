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
    
public:
    
    void load(UserInput userInput);
    
    void read(InSequences& inSequence);
    
};

#endif /* GFASTATS_INPUT_H */
