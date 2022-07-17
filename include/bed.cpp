#include <stdlib.h>
#include <vector>
#include <string>

#include "bed.h"

void BedCoordinates::pushCoordinates(std::string h, unsigned int b, unsigned int e) { // reading coordinates
    
    seqHeaders.push_back(h);
    cBegin.push_back(b);
    cEnd.push_back(e);
    
}

bool BedCoordinates::empty() {
    
    return (seqHeaders.size()==0) ? true : false; // check if no coordinates are present
    
}

unsigned int BedCoordinates::size() {
    
    return seqHeaders.size(); // check if no coordinates are present
    
}

std::vector<std::string> BedCoordinates::getSeqHeaders() { // get all the headers
    
    return seqHeaders;
    
}

std::string BedCoordinates::getSeqHeader(unsigned int pos) { // get a specific header
    
    return seqHeaders[pos];
    
}

unsigned int BedCoordinates::getcBegin(unsigned int pos) { // get a specific start coordinate
    
    return cBegin[pos];
    
}

unsigned int BedCoordinates::getcEnd(unsigned int pos) { // get a specific end coordinate
    
    return cEnd[pos];
    
}
