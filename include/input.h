#ifndef INPUT_H
#define INPUT_H

struct UserInputGfastats : UserInput {
    
    std::vector<std::string> outFiles; // output files
    int seqReport_flag = 0;
    int outSequence_flag = 0;
    int nstarReport_flag = 0;
    int outSize_flag = 0;
    int outCoord_flag = 0;
    int outFile_flag = 0;
    int outBubbles_flag = 0;
    int cmd_flag = 0;
    int rmGaps_flag = 0;
    int extractContigs_flag = 0;
    int terminalOvlLen = 0;

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
    
public:
    
    void load(UserInputGfastats userInput);
    
    void read(InSequences& inSequence);
    
};

#endif /* INPUT_H */
