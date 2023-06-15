#ifndef INPUT_H
#define INPUT_H

struct UserInputGfastats : UserInput {
    
    std::vector<std::string> outFiles; // output files

};

class Input {
    
    UserInputGfastats userInput;
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    std::shared_ptr<std::istream> stream;
    
    std::vector<Instruction> instructions;
    
    unsigned int begin = 0, end = 0;
    
public:
    
    void load(UserInputGfastats userInput);
    
    void read(InSequences& inSequence);
    
};

#endif /* INPUT_H */
