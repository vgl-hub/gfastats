#include <stdlib.h>
#include <string>
#include <vector>
#include "gfastats-struct.h"
#include "gfastats-functions.h"
#include "gfa-lines.h"
    
void InSegment::setSeqHeader(std::string* h) {
    seqHeader = *h;
}

void InSegment::setSeqComment(std::string c) {
    seqComment = c;
}

void InSegment::setInSequence(std::string* s) {
    inSequence = *s;
}

void InSegment::setInSequenceQuality(std::string* q) {
    inSequenceQuality = *q;
}

void InSegment::setSeqTags(std::vector<Tag>* t) {
    tags = *t;
}

void InSegment::setuId(unsigned int i) { // absolute id
    uId = i;
}

void InSegment::setiId(unsigned int i) { // temporary id, internal to scaffold
    iId = i;
}

void InSegment::setSeqPos(unsigned int i) { // temporary id, internal to scaffold
    seqPos = i;
}

std::string InSegment::getSeqHeader() {
    return seqHeader;
}

std::string InSegment::getSeqComment() {
    return seqComment;
}

std::vector<Tag> InSegment::getTags() {
    return tags;
}

std::string InSegment::getInSequence(unsigned int start, unsigned int end) {
    
    if (inSequence == "") {
        
        return "*";
        
    }else{
    
        return start != 0 || end != 0 ? inSequence.substr(start-1, end-start+1) : inSequence;
        
    }
    
}

std::string InSegment::getInSequenceQuality(unsigned int start, unsigned int end) {
    
    return start != 0 || end != 0 ? inSequenceQuality.substr(start-1, end-start+1) : inSequenceQuality;
    
}

unsigned int InSegment::getSeqPos() {
    
    return seqPos;

}

unsigned long long int InSegment::getSegmentLen(unsigned long long int start, unsigned long long int end) {
    
    if (inSequence == "") {
        
        return lowerCount;
        
    }else{
    
        return start != 0 || end != 0 ? end-start+1 : A + C + G + T; // need to sum long long int to prevent size() overflow
        
    }
    
}

unsigned int InSegment::getuId() { // absolute id
    
    return uId;
}

unsigned int InSegment::getiId() { // temporary id, internal to scaffold
    
    return iId;
}

void InSegment::setACGT(unsigned long long int* a, unsigned long long int* c, unsigned long long int* g, unsigned long long int* t) {
    
    A = *a;
    C = *c;
    G = *g;
    T = *t;
    
}

void InSegment::setLowerCount(unsigned long long int* C) {
    
    lowerCount = *C;
    
}

unsigned long long int InSegment::getA() {
    
    return A;
}

unsigned long long int InSegment::getC() {
    
    return C;
}

unsigned long long int InSegment::getG() {
    
    return G;
}

unsigned long long int InSegment::getT() {
    
    return T;
}

unsigned int InSegment::getLowerCount(unsigned long long int start, unsigned long long int end) {
    
    if (start == 0 || end == 0) {
        
        return lowerCount;
        
    }else{
        
        unsigned long long int lowerCountSubset = 0;
        
        for (char base : inSequence) { // need to fix this loop
            
            if (islower(base)) {
                
                ++lowerCountSubset;
                
            }
            
        }
        
        return lowerCountSubset;
        
    }

}

double InSegment::computeGCcontent() {
    
    double GCcontent = (double) (G + C) / (G + C + A + T) * 100;
    
    return GCcontent;
}

bool InSegment::trimSegment(unsigned int start, unsigned int end) {
    
    for(char& base : inSequence.substr(start, end-start)) {
        
        switch (base) {
            case 'A':
            case 'a':{
                
                A--;
                break;
                
            }
            case 'C':
            case 'c':{
                
                C--;
                break;
                
            }
            case 'G':
            case 'g': {
                
                G--;
                break;
                
            }
            case 'T':
            case 't': {
                
                T--;
                break;
                
            }
                
        }
        
    }
    
    inSequence.erase(start, end-start);
    
    if (inSequenceQuality.size()>0) {
    
        inSequenceQuality.erase(start, end-start);
    
    }
    
    return true;
}

bool InSegment::rvcpSegment() {

    inSequence = revCom(inSequence);

    return true;
    
}

bool InSegment::invertSegment() {

    inSequence = rev(inSequence);
    inSequenceQuality = rev(inSequenceQuality);

    return true;
    
}

void InGap::newGap(unsigned int uId, unsigned int sId1, unsigned int sId2, const char& sId1or, const char& sId2or, unsigned int& dist, std::string gHeader, std::vector<Tag> tags) {
    
    this->gHeader = gHeader;
    this->uId = uId;
    this->sId1 = sId1;
    this->sId2 = sId2;
    this->sId1Or = sId1or;
    this->sId2Or = sId2or;
    this->dist = dist;
    this->tags = tags;
    
}

void InGap::setuId(unsigned int i) { // absolute id
    uId = i;
}

void InGap::setiId(unsigned int i) { // temporary id, internal to scaffold
    iId = i;
}

void InGap::setsId1(unsigned int i) {
    sId1 = i;
}

void InGap::setsId2(unsigned int i) {
    sId2 = i;
}

void InGap::setDist(unsigned int i) {
    dist = i;
}

std::string InGap::getgHeader() {
    
    return gHeader;
    
}

unsigned int InGap::getuId() {
    
    return uId;
    
}

unsigned int InGap::getsId1() {
    
    return sId1;
    
}

char InGap::getsId1Or() {
    
    return sId1Or;
    
}

unsigned int InGap::getsId2() {
    
    return sId2;
    
}

char InGap::getsId2Or() {
    
    return sId2Or;
    
}

unsigned int InGap::getDist(unsigned int start, unsigned int end) {
    
    return start != 0 || end != 0 ? end-start+1 : dist;
    
}

std::vector<Tag> InGap::getTags() {
    
    return tags;
    
}

void InEdge::newEdge(unsigned int eUId, unsigned int sId1, unsigned int sId2, const char& sId1Or, const char& sId2Or, std::string cigar, std::string eHeader, std::vector<Tag> tags) {
    
    this->eUId = eUId;
    this->sId1 = sId1;
    this->sId2 = sId2;
    this->sId1Or = sId1Or;
    this->sId2Or = sId2Or;
    this->cigar = cigar;
    this->eHeader = eHeader;
    this->tags = tags;
    
}

bool InEdge::operator==(const InEdge& e) const {
    return sId1 == e.sId1 && sId2 == e.sId2 && sId1Or == e.sId1Or && e.sId2Or == sId2Or;
}

void InEdge::seteUId(unsigned int i) { // absolute id
    eUId = i;
}

void InEdge::seteId(unsigned int i) { // temporary id, internal to scaffold
    eId = i;
}

void InEdge::setsId1(unsigned int i) {
    sId1 = i;
}

void InEdge::setsId2(unsigned int i) {
    sId2 = i;
}

void InEdge::setSeqTags(std::vector<Tag>* t) {
    tags = *t;
}

std::string InEdge::getCigar() {
    
    return cigar;
    
}

unsigned int InEdge::geteUId() {
    
    return eUId;
    
}

unsigned int InEdge::geteId() {
    
    return eId;
    
}

unsigned int InEdge::getsId1() {
    
    return sId1;
    
}

char InEdge::getsId1Or() {
    
    return sId1Or;
    
}

unsigned int InEdge::getsId2() {
    
    return sId2;
    
}

char InEdge::getsId2Or() {
    
    return sId2Or;

}

std::vector<Tag> InEdge::getTags() {
    return tags;
}
    
void InPath::newPath(unsigned int pUid, std::string h, std::string c, unsigned int seqpos) {
    
    pHeader = h;
    pComment = c;
    pathComponents.clear();
    pUId = pUid;
    seqPos = seqpos;

}

void InPath::setpUId(unsigned int pUid) {
    
    pUId = pUid;

}

void InPath::setHeader(std::string pheader) {
    
    pHeader = pheader;

}

void InPath::setComment(std::string c) {
    pComment = c;
}

void InPath::add(PathType type, unsigned int UId, char sign, unsigned long long int start, unsigned long long int end) {
    
    pathComponents.push_back({type, UId, sign, start, end});
    
}

void InPath::append(std::vector<PathComponent> components) {

    pathComponents.insert(std::end(pathComponents), std::begin(components), std::end(components));
    
}

void InPath::clearPath() {
    
    pathComponents.clear();
    
}

void InPath::setComponents(std::vector<PathComponent> newComponents) {

    pathComponents = newComponents;
    
}

std::vector<PathComponent> InPath::getComponents() {
    
    return pathComponents;
    
}

std::vector<PathComponent>* InPath::getComponentsByRef() {
    
    return &pathComponents;
    
}

unsigned int InPath::getpUId() {
    
    return pUId;
    
}

std::string InPath::getHeader() {
    
    return pHeader;
    
}

std::string InPath::getComment() {
    
    return pComment;

}

unsigned int InPath::getSeqPos() {
    
    return seqPos;

}

unsigned int InPath::getContigN() {
    
    return contigN;
    
}

unsigned long long int InPath::getLen() {
    
    return length;
    
}

unsigned long long int InPath::getA() {
    
    return A;
    
}

unsigned long long int InPath::getC() {
    
    return C;
    
}

unsigned long long int InPath::getG() {
    
    return G;
    
}

unsigned long long int InPath::getT() {
    
    return T;
    
}

unsigned long long int InPath::getSegmentLen() {
    
    return segmentLength;
    
}

unsigned long long int InPath::getLowerCount() {
    
    return lowerCount;
    
}

void InPath::revCom() {
    
    revComPathComponents(pathComponents);

}

void InPath::increaseContigN() {
    
    contigN++;

}

void InPath::increaseGapN() {
    
    contigN++;

}

void InPath::increaseLen(unsigned long long int n) {
    
    length += n;

}

void InPath::increaseSegmentLen(unsigned long long int n) {
    
    segmentLength += n;

}

void InPath::increaseLowerCount(unsigned long long int n) {
    
    lowerCount += n;

}

void InPath::increaseA(unsigned long long int n) {
    
    A += n;

}

void InPath::increaseC(unsigned long long int n) {
    
    C += n;

}

void InPath::increaseG(unsigned long long int n) {
    
    G += n;

}

void InPath::increaseT(unsigned long long int n) {
    
    T += n;

}
