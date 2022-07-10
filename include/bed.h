#ifndef BED_H
#define BED_H

class BedCoordinates { // generic representation of bed coordinates
private:
    std::vector<std::string> seqHeaders;
    std::vector<unsigned int> cBegin;
    std::vector<unsigned int> cEnd;
    
public:
    
    void pushCoordinates(std::string h, unsigned int b = 0, unsigned int e = 0); // reading coordinates
    
    bool empty(); // check if no coordinates are present
    
    unsigned int size();
    
    std::vector<std::string> getSeqHeaders(); // get all the headers
    
    std::string getSeqHeader(unsigned int pos); // get a specific header
    
    unsigned int getcBegin(unsigned int pos); // get a specific start coordinate
    
    unsigned int getcEnd(unsigned int pos); // get a specific end coordinate
    
};

#endif /* BED_H */
