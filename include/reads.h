#ifndef READS_H
#define READS_H

class InReads {
    
    std::vector<Log> logs;
    
    UserInput userInput;
    std::vector<InSegment*> inReads;
    
    //intermediates
    std::string h;
    char* c;
    
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::vector<unsigned long long int> readNstars    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> readLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
public:
    
    ~InReads();
    
    void load(UserInput userInput);
    
    void traverseInReads(Sequences* sequence);
    
    InSegment* traverseInRead(Log* threadLog, Sequence* sequence, unsigned int seqPos);
    
    void appendReads(Sequences* sequence);
    
    unsigned long long int getTotReadLen();
    
    double computeAvgReadLen();
    
    unsigned long long int getReadN50();
    
    void report();
    
    void evalNstars();
    
};

#endif /* READS_H */
