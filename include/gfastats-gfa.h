//
//  gfastats-gfa.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFASTATS_GFA_H
#define GFASTATS_GFA_H

//classes
class BedCoordinates { // generic representation of bed coordinates
private:
    std::vector<std::string> seqHeaders;
    std::vector<unsigned int> cBegin;
    std::vector<unsigned int> cEnd;
    
public:
    
    void pushCoordinates(std::string h, unsigned int b = 0, unsigned int e = 0) { // reading coordinates
        
        seqHeaders.push_back(h);
        cBegin.push_back(b);
        cEnd.push_back(e);
        
    }
    
    bool empty() {
        
        return (seqHeaders.size()==0) ? true : false; // check if no coordinates are present
        
    }
    
    unsigned int size() {
        
        return seqHeaders.size(); // check if no coordinates are present
        
    }
    
    std::vector<std::string> getSeqHeaders() { // get all the headers
        
        return seqHeaders;
        
    }
    
    std::string getSeqHeader(unsigned int pos) { // get a specific header
        
        return seqHeaders[pos];
        
    }
    
    unsigned int getcBegin(unsigned int pos) { // get a specific start coordinate
        
        return cBegin[pos];
        
    }
    
    unsigned int getcEnd(unsigned int pos) { // get a specific end coordinate
        
        return cEnd[pos];
        
    }
    
};

class InSegment { // DNA sequence with no gaps
private:
    std::string seqHeader;
    std::string seqComment;
    std::string inSequence;
    std::string inSequenceQuality;
    unsigned int A = 0, C = 0, G = 0, T = 0, lowerCount = 0, uId = 0, iId = 0;
    
    friend class SAK;
    friend class InSequences;
    
public:
    
    void setSeqHeader(std::string* h) {
        seqHeader = *h;
    }
    
    void setSeqComment(std::string c) {
        seqComment = c;
    }
    
    void setInSequence(std::string* s) {
        inSequence = *s;
    }
    
    void setInSequenceQuality(std::string* q) {
        inSequenceQuality = *q;
    }

    void setuId(unsigned int i) { // absolute id
        uId = i;
    }
    
    void setiId(unsigned int i) { // temporary id, internal to scaffold
        iId = i;
    }
    
    std::string getSeqHeader() {
        return seqHeader;
    }
    
    std::string getSeqComment() {
        return seqComment;
    }
    
    std::string getInSequence(unsigned int start = 0, unsigned int end = 0) {
        
        return start != 0 || end != 0 ? inSequence.substr(start-1, end-start+1) : inSequence;
        
    }
    
    std::string getInSequenceQuality(unsigned int start = 0, unsigned int end = 0) {
        
        return start != 0 || end != 0 ? inSequenceQuality.substr(start-1, end-start+1) : inSequenceQuality;
        
    }
    
    unsigned int getSegmentLen() {
        return inSequence.size();
    }
    
    unsigned int getuId() { // absolute id
        
        return uId;
    }
    
    unsigned int getiId() { // temporary id, internal to scaffold
        
        return iId;
    }
    
    void setACGT(unsigned int* a, unsigned int* c, unsigned int* g, unsigned int* t) {
        
        A = *a;
        C = *c;
        G = *g;
        T = *t;
        
    }
    
    void setLowerCount(unsigned int* C) {
        
        lowerCount = *C;
        
    }
    
    unsigned int getA() {
        
        return A;
    }
    
    unsigned int getC() {
        
        return C;
    }
    
    unsigned int getG() {
        
        return G;
    }
    
    unsigned int getT() {
        
        return T;
    }
    
    unsigned int getLowerCount(unsigned int start = 0, unsigned int end = 0) {
        
        if (start == 0 || end == 0) {
            
            return lowerCount;
            
        }else{
            
            long unsigned int lowerCountSubset = 0;
            
            for (char base : inSequence) {
                
                if (islower(base)) {
                    
                    ++lowerCountSubset;
                    
                }
                
            }
            
            return lowerCountSubset;
            
        }

    }
    
    unsigned int getSegmentLength(unsigned int start = 0, unsigned int end = 0) {
        
        return start != 0 || end != 0 ? end-start+1 : inSequence.size();
        
    }
    
    double computeGCcontent() {
        
        double GCcontent = (double) (G + C) / (G + C + A + T) * 100;
        
        return GCcontent;
    }
    
    bool trimSegment(unsigned int start, unsigned int end) {
        
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
    
    bool rvcpSegment() {
    
        inSequence = revCom(inSequence);
    
        return true;
        
    }
    
    bool invertSegment() {
    
        inSequence = rev(inSequence);
        inSequenceQuality = rev(inSequenceQuality);
    
        return true;
        
    }
    
};

class InGap {
private:
//    unsigned long long int lineN; // useful if we wish to sort as is the original input
    std::string gHeader;
    char sId1Or, sId2Or;
    unsigned int uId, iId, sId1, sId2, dist;
    
    friend class SAK;
    friend class InSequences;
    
public:
    void newGap(unsigned int uid, unsigned int sid1, unsigned int sid2, const char& sid1or, const char& sid2or, unsigned int& d, std::string gheader = "") {
        
        gHeader = gheader;
        uId = uid;
        sId1 = sid1;
        sId2 = sid2;
        sId1Or = sid1or;
        sId2Or = sid2or;
        dist = d;
        
    }

    void setuId(unsigned int i) { // absolute id
        uId = i;
    }
    
    void setiId(unsigned int i) { // temporary id, internal to scaffold
        iId = i;
    }

    void setsId1(unsigned int i) {
        sId1 = i;
    }
    
    void setsId2(unsigned int i) {
        sId2 = i;
    }
    
    void setDist(unsigned int i) {
        dist = i;
    }
    
    std::string getgHeader() {
        
        return gHeader;
        
    }
    
    unsigned int getuId() {
        
        return uId;
        
    }
    
    unsigned int getsId1() {
        
        return sId1;
        
    }
    
    char getsId1Or() {
        
        return sId1Or;
        
    }
    
    unsigned int getsId2() {
        
        return sId2;
        
    }
    
    char getsId2Or() {
        
        return sId2Or;
        
    }
    
    unsigned int getDist(unsigned int start = 0, unsigned int end = 0) {
        
        return start != 0 || end != 0 ? end-start+1 : dist;
        
    }
    
};
class InEdge {
    private:
//    unsigned long long int lineN; // useful if we wish to sort as is the original input
    std::string cigar, eHeader;
    char sId1Or, sId2Or;
    unsigned int eUId, eId, sId1, sId2;
    
    friend class SAK;
    friend class InSequences;
    
public:
    void newEdge(unsigned int eUid, unsigned int sid1, unsigned int sid2, const char& sid1or, const char& sid2or, std::string c = "", std::string h = "") {
        
        eUId = eUid;
        sId1 = sid1;
        sId2 = sid2;
        sId1Or = sid1or;
        sId2Or = sid2or;
        cigar = c;
        eHeader = h;
        
    }

    void seteUId(unsigned int i) { // absolute id
        eUId = i;
    }
    
    void seteId(unsigned int i) { // temporary id, internal to scaffold
        eId = i;
    }

    void setsId1(unsigned int i) {
        sId1 = i;
    }
    
    void setsId2(unsigned int i) {
        sId2 = i;
    }
    
    std::string getCigar() {
        
        return cigar;
        
    }
    
    unsigned int geteUId() {
        
        return eUId;
        
    }

    unsigned int geteId() {
        
        return eId;
        
    }
    
    unsigned int getsId1() {
        
        return sId1;
        
    }
    
    char getsId1Or() {
        
        return sId1Or;
        
    }
    
    unsigned int getsId2() {
        
        return sId2;
        
    }
    
    char getsId2Or() {
        
        return sId2Or;
    
    }
};

class InPath {
    
private:
//    unsigned long long int lineN; // useful if we wish to sort as is the original input
    std::string pHeader, pComment;
    std::vector<PathTuple> pathComponents;
    unsigned int pUId;
    
    long unsigned int length = 0, lowerCount = 0, A = 0, C = 0, G = 0, T = 0;
    
    friend class SAK;
    friend class InSequences;

public:
    
    void newPath(unsigned int pUid, std::string h, std::string c = "") {
        
        pHeader = h;
        pComment = c;
        pathComponents.clear();
        pUId = pUid;
        
        verbose("Processed sequence: " + pHeader + " (uId: " + std::to_string(pUId) + ")");
    
    }

    void setpUId(unsigned int pUid) {
        
        pUId = pUid;
    
    }
    
    void setHeader(std::string pheader) {
        
        pHeader = pheader;
    
    }
    
    void setComment(std::string c) {
        pComment = c;
    }
    
    void add(char type, unsigned int UId, char sign = '+', unsigned int start = 0, unsigned int end = 0) {
        
        pathComponents.push_back(std::make_tuple(type, UId, sign, start, end));
        
    }
    
    void append(std::vector<PathTuple> components) {
    
        pathComponents.insert(std::end(pathComponents), std::begin(components), std::end(components));
        
    }
    
    void clearPath() {
        
        pathComponents.clear();
        
    }
    
    void setComponents(std::vector<PathTuple> newComponents) {

        pathComponents = newComponents;
        
    }
    
    std::vector<PathTuple> getComponents() {
        
        return pathComponents;
        
    }
    
    std::vector<PathTuple>* getComponentsByRef() {
        
        return &pathComponents;
        
    }

    unsigned int getpUId() {
        
        return pUId;
        
    }
    
    std::string getHeader() {
        
        return pHeader;
        
    }
    
    std::string getComment() {
        
        return pComment;
    
    }
    
    unsigned int getLen() {
        
        return length;
        
    }
    
    unsigned int getA() {
        
        return A;
        
    }
    
    unsigned int getC() {
        
        return C;
        
    }
    
    unsigned int getG() {
        
        return G;
        
    }
    
    unsigned int getT() {
        
        return T;
        
    }
    
    unsigned int getLowerCount() {
        
        return lowerCount;
        
    }
    
    void revCom() {
        
        revComPathComponents(pathComponents);
    
    }
    
    void increaseLen(unsigned int n) {
        
        length += n;
    
    }
    
    void increaseLowerCount(unsigned int n) {
        
        lowerCount += n;
    
    }
    
    void increaseA(unsigned int n) {
        
        A += n;
    
    }
    
    void increaseC(unsigned int n) {
        
        C += n;
    
    }
    
    void increaseG(unsigned int n) {
        
        G += n;
    
    }
    
    void increaseT(unsigned int n) {
        
        T += n;
    
    }
    
};

class InSequences { //collection of InSegment and inGap objects and their summary statistics
    
private:
    
    //gfa variables
    std::vector<InSegment> inSegments;
    std::vector<InGap> inGaps;
    std::vector<InEdge> inEdges;
    std::vector<InPath> inPaths;
    std::vector<std::vector<Tuple>> adjListFW;
    std::vector<std::vector<Tuple>> adjListBW;
    std::vector<std::vector<EdgeTuple>> adjEdgeListFW;
    phmap::flat_hash_map<std::string, unsigned int> headersToIds;
    phmap::flat_hash_map<unsigned int, std::string> idsToHeaders;
    phmap::flat_hash_map<int, bool> visited, deleted;
    bool backward = false;
    
    std::vector<unsigned int> scaffLens;
    std::vector<unsigned int> contigLens;
    std::vector<unsigned int> gapLens;
    
    std::vector<unsigned int> scaffNstars   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffLstars   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffNGstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffLGstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<unsigned int> contigNstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigLstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigNGstars {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigLGstars {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<unsigned int> gapNstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> gapLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    double scaffAuN = 0, scaffAuNG = 0, contigAuN = 0, contigAuNG = 0, gapAuN = 0;
    
    InSegment inSegment;
    InGap gap;
    InEdge edge;
    InPath path;
    
    unsigned long long int totScaffLen = 0, totSegmentLen = 0;
    
    unsigned int
    scaffN = 0,
    pathN = 0,
    totGapLen = 0,
    uId = 0; // unique numeric identifier for each feature
    
    unsigned long int totA = 0;
    unsigned long int totC = 0;
    unsigned long int totG = 0;
    unsigned long int totT = 0;
    unsigned long int totLowerCount = 0;

    //connectivity
    unsigned int deadEnds = 0;
    unsigned int disconnectedComponents = 0;
    unsigned int lengthDisconnectedComponents = 0;
    
    std::vector<Bubble> bubbles;
    
    friend class SAK;
    
public:
    
    void addSegment(unsigned int uId, unsigned int iId, std::string seqHeader, std::string* seqComment, std::string* sequence, unsigned int* A, unsigned int* C, unsigned int* G, unsigned int* T, unsigned int* lowerCount, std::string* sequenceQuality = NULL) {
        
        // operations on the segment
        
        inSegment.setiId(iId); // set temporary sId internal to scaffold
        
        inSegment.setuId(uId); // set absolute id
        
        inSegment.setSeqHeader(&seqHeader);
        
        if (seqComment != NULL) {
            
            inSegment.setSeqComment(*seqComment);
            
        }
        
        verbose("Processing segment: " + seqHeader + " (uId: " + std::to_string(uId) + ", iId: " + std::to_string(iId) + ")");
        
        inSegment.setInSequence(sequence);
        
        verbose("Segment sequence set");
        
        if (sequenceQuality != NULL) {
            
            inSegment.setInSequenceQuality(sequenceQuality);
            
            verbose("Segment sequence quality set");
            
        }
        
        inSegment.setACGT(A, C, G, T);
        
        verbose("Increased ACGT counts");
        
        totA += *A;
        totC += *C;
        totG += *G;
        totT += *T;
        
        verbose("Increased total ACGT counts");
        
        inSegment.setLowerCount(lowerCount);
        
        verbose("Increased count of lower bases");
        
        totLowerCount += *lowerCount;

        verbose("Increased total count of lower bases");
        
        inSegments.push_back(inSegment); // adding segment to segment set
        
        // operations of the segment set
        
        verbose("Segment added to sequence vector");
        
        unsigned int seqSize = sequence->size();
        
        contigLens.push_back(seqSize);
        
        verbose("Recorded length of sequence");
        
        changeTotSegmentLen(seqSize);
        
        verbose("Increased total segment length");
        
    }
    
    void pushbackGap(std::string* seqHeader, unsigned int* iId, unsigned int* uId, unsigned int pos, unsigned int* dist, char sign, unsigned int seqLen, int n) {
        
        verbose("Processing gap " + *seqHeader+"."+std::to_string(*iId) + " (uId: " + std::to_string(*uId) + ", iId: " + std::to_string(*iId) + ")");
        
        gap.newGap(*uId, (pos - *dist == 0) ? *uId+1 : *uId-1, (pos - seqLen == 0) ? *uId+n : *uId+1, '+', sign, *dist, *seqHeader+"."+std::to_string(*iId));
        
        addGap(gap);
        
        insertHash1(*seqHeader+"."+std::to_string(*iId), *uId); // header to hash table
        insertHash2(*uId, *seqHeader+"."+std::to_string(*iId)); // uID to hash table
        
        path.add('G', *uId, '0');
        
        *dist=0;
        
        (*iId)++; // number of gaps in the current scaffold
        (*uId)++; // unique numeric identifier
        
    }
    
    void pushbackSegment(std::string* seqHeader, std::string* seqComment, std::string* sequence, unsigned int* iId, unsigned int* uId, unsigned int* A, unsigned int* C, unsigned int* G, unsigned int* T, unsigned int* lowerCount, unsigned int sStart, unsigned int sEnd, std::string* sequenceQuality = NULL) {
         
        std::string sequenceSubSeq, sequenceQualitySubSeq;
        
        sequenceSubSeq = sequence->substr(sStart, sEnd + 1 - sStart);
        
        if (sequenceQuality != NULL) {
            
            sequenceQualitySubSeq = sequenceQuality->substr(sStart, sEnd + 1 - sStart);
            
        }
        
        addSegment(*uId, *iId, *seqHeader+"."+std::to_string(*iId), seqComment, &sequenceSubSeq, A, C, G, T, lowerCount, &sequenceQualitySubSeq);
        
        insertHash1(*seqHeader+"."+std::to_string(*iId), *uId); // header to hash table
        insertHash2(*uId, *seqHeader+"."+std::to_string(*iId)); // uID to hash table
        
        path.add('S', *uId, '+');
        
        *A = 0, *C = 0, *G = 0, *T = 0, *lowerCount = 0;
        
        (*iId)++; // number of gaps in the current scaffold
        (*uId)++; // unique numeric identifier
        
    }

    void homopolymerCompress(std::string *sequence, std::vector<unsigned int> &compressionIndices, std::vector<unsigned int> &compressionLengths, unsigned int cutoff) {
        unsigned int index=0, length, new_length=0;

        auto lambda = [&length, &index, &compressionIndices, &compressionLengths, &sequence, &new_length, &cutoff](int i){
            length = i-index;
            if(length > cutoff) {
                compressionIndices.push_back(new_length);
                compressionLengths.push_back(length);
            }
            int num = length > cutoff ? 1 : length;
            memset(&((*sequence)[new_length]), (*sequence)[index], num);
            new_length += num;
        };

        for(unsigned int i=1; i<sequence->length(); ++i) {
            if((*sequence)[i] == (*sequence)[index]) continue;
            lambda(i);
            index = i;
        }
        lambda(sequence->length());

        sequence->resize(new_length);
    }

    void homopolymerDecompress(std::string *sequence, const std::vector<unsigned int> &compressionIndices, const std::vector<unsigned int> &compressionLengths) {
        std::string ret="";
        ret.reserve(sequence->length()*2); // random guess for final sequence length
        for(unsigned int i=0, ci=0, len; i<sequence->length(); ++i) {
            len = ((ci<compressionIndices.size() && i==compressionIndices[ci]) ? compressionLengths[ci++] : 1);
            ret += std::string(len, (*sequence)[i]);
        }
        ret.shrink_to_fit();
        *sequence = ret;
    }
    
    void traverseInSequence(std::string* pHeader, std::string* seqComment, std::string* sequence, std::string* sequenceQuality = NULL) { // traverse the sequence to split at gaps and measure sequence properties
        
        std::vector<unsigned int> compressionIndices, compressionLengths;
        if(hc_flag) {
            homopolymerCompress(sequence, compressionIndices, compressionLengths, hc_cutoff);
        }

        unsigned int pos = 0, // current position in sequence
        hc_index=0, // used with homopolymer compression
        A = 0, C = 0, G = 0, T = 0,
        lowerCount = 0,
        dist = 0, // gap size
        iId = 1, // scaffold feature internal identifier
        sStart = 0, sEnd = 0; // segment start and end
        char sign = '+';
        bool wasN = false;
        
        path.clearPath();
        
        phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find (*pHeader); // get the headers to uIds table to look for the header
        
        if (got == headersToIds.end()) { // this is the first time we see this path name
            
            insertHash1(*pHeader, uId); // header to hash table
            insertHash2(uId, *pHeader); // uId to hash table
            
        }else{
            
            fprintf(stderr, "Error: path name already exists (%s). Terminating.\n", pHeader->c_str()); exit(1);
            
        }
        
        path.newPath(uId, *pHeader);
        
        uId++;
        
        if (seqComment != NULL) {
            
            path.setComment(*seqComment);
            
        }
        
        unsigned int seqLen = sequence->length()-1;
        for (char &base : *sequence) {
            unsigned int count = hc_flag && hc_index < compressionIndices.size() && pos == compressionIndices[hc_index] ? compressionLengths[hc_index++] : 1;
            
            if (islower(base)) {
                
                lowerCount+=count;
                
            }
            
            switch (base) {
                    
                case 'N':
                case 'n':
                case 'X':
                case 'x': {
                    
                    dist+=count;
                    
                    if (!wasN && pos>0) { // gap start and gap not at the start of the sequence
                            
                        sEnd = pos - 1;
                        pushbackSegment(pHeader, seqComment, sequence, &iId, &uId, &A, &C, &G, &T, &lowerCount, sStart, sEnd, sequenceQuality);
                        
                    }
                    
                    if(pos == seqLen) { // end of scaffold, terminal gap
                        
                        sign = '-';
                        
                        pushbackGap(pHeader, &iId, &uId, pos, &dist, sign, seqLen, -1);
                        
                    }
                    
                    wasN = true;
                    
                    break;
                }
                default: {
                    
                    switch (base) {
                        case 'A':
                        case 'a':{
                            
                            A+=count;
                            break;
                            
                        }
                        case 'C':
                        case 'c':{
                            
                            C+=count;
                            break;
                            
                        }
                        case 'G':
                        case 'g': {
                            
                            G+=count;
                            break;
                            
                        }
                        case 'T':
                        case 't': {
                            
                            T+=count;
                            break;
                            
                        }
                            
                    }
                    
                    if (wasN) { // internal gap end
                        
                        sStart = pos;
                        pushbackGap(pHeader, &iId, &uId, pos, &dist, sign, seqLen, 1);
                        
                    }
                    
                    if (pos == seqLen) {
                        
                        sEnd = pos;
                        pushbackSegment(pHeader, seqComment, sequence, &iId, &uId, &A, &C, &G, &T, &lowerCount, sStart, sEnd, sequenceQuality);
                        
                    }
                    
                    wasN = false;
                    
                }
                    
            }
            
            pos++;
            
        }
        
        inPaths.push_back(path);
        
        verbose("Added fasta sequence as path");
        
        recordScaffLen(seqLen+1);
        
        verbose("Recorded length of sequence");
        
    }
    
    void traverseInSegment(std::string* seqHeader, std::string* seqComment, std::string* sequence, std::string* sequenceQuality = NULL) { // traverse the sequence to split at gaps and measure sequence properties
        
        unsigned int A = 0, C = 0, G = 0, T = 0, sUId = 0, lowerCount = 0;
        
        for (char &base : *sequence) {
            
            if (islower(base)) {
                
                lowerCount++;
                
            }
                    
            switch (base) {
                case 'A':
                case 'a':{
                    
                    A++;
                    break;
                    
                }
                case 'C':
                case 'c':{
                    
                    C++;
                    break;
                    
                }
                case 'G':
                case 'g': {
                    
                    G++;
                    break;
                    
                }
                case 'T':
                case 't': {
                    
                    T++;
                    break;
                    
                }
                    
            }
                
        }
        
        phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find (*seqHeader); // get the headers to uIds table to look for the header
        
        if (got == headersToIds.end()) { // this is the first time we see this segment
            
            insertHash1(*seqHeader, uId); // header to hash table
            insertHash2(uId, *seqHeader); // header to hash table
            
            sUId = uId;
            
            uId++;
            
        }else{
            
            sUId = got->second;
            
        }
                
        addSegment(sUId, 0, *seqHeader, seqComment, sequence, &A, &C, &G, &T, &lowerCount, sequenceQuality);
        
        recordScaffLen(sequence->size());
        
        verbose("Recorded length of sequence");
        
    }
    
    void appendSequence(std::string* pHeader, std::string* pComment, std::string* sequence, std::string* sequenceQuality = NULL) { // method to append a new sequence from a fasta
        
        verbose("Header, comment, sequence and (optionally) quality read");
        
        if(verbose_flag) {std::cerr<<"\n";};
        
        traverseInSequence(pHeader, pComment, sequence, sequenceQuality);
        
        verbose("Sequence traversed");
        
        scaffN++;
        
        verbose("Increased total scaffold N");
        
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    void appendSegment(std::string* seqHeader, std::string* seqComment, std::string* sequence, std::string* sequenceQuality = NULL) { // method to append a new segment from a gfa
        
        verbose("Header, comment, sequence and (optionally) quality read");
        
        if(verbose_flag) {std::cerr<<"\n";};
        
        traverseInSegment(seqHeader, seqComment, sequence, sequenceQuality);
        
        verbose("Segment traversed");
        
        scaffN++;
        
        verbose("Increased total scaffold N");
        
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    unsigned int getSegmentN() {
        
        return inSegments.size();
    }
    
    InSegment getInSegment(unsigned int sId) {
        
        auto inSegment = find_if(inSegments.begin(), inSegments.end(), [sId](InSegment& obj) {return obj.getuId() == sId;}); // given a uId, find it in nodes
        
        if (inSegment == inSegments.end()) {
        
            fprintf(stderr, "Error: could not find segment (sId: %i).\n", sId); exit(1);
            
        }
            
        return *inSegment;
        
    }
    
    std::vector<InSegment>* getInSegments() {
        
        return &inSegments;
        
    }
    
    std::vector<InGap>* getInGaps() {
        
        return &inGaps;
        
    }
    
    InPath getInPath(unsigned int pId) {
        
        auto inPath = find_if(inPaths.begin(), inPaths.end(), [pId](InPath& obj) {return obj.getpUId() == pId;}); // given a uId, find it in nodes
        
        if (inPath == inPaths.end()) {
        
            fprintf(stderr, "Error: could not find path (pId: %i).\n", pId); exit(1);
            
        }
            
        return *inPath;
        
    }
    
    std::vector<InPath> getInPaths() {
        
        return inPaths;
        
    }
    
    unsigned long long int getTotScaffLen() {
        
        return totScaffLen;
        
    }
    
    void changeTotSegmentLen(int segmentLen) {
        
        totSegmentLen += segmentLen;
        
    }
    
    unsigned long long int getTotSegmentLen() {
        
        return totSegmentLen;
        
    }
    
    void changeTotGapLen(unsigned int gapLen) {
        
        totGapLen += gapLen;
        
    }
    
    unsigned int getGapN() {
        
        return gapLens.size();
        
    }

    unsigned int getEdgeN() {
        
        return inEdges.size();
        
    }
    
    unsigned int getPathN() {
        
        return inPaths.size();
        
    }
    
    void recordScaffLen(unsigned int seqLen) {
        
        scaffLens.push_back(seqLen);
        
    }
    
    void recordGapLen(unsigned int gapLen) {
        
        gapLens.push_back(gapLen);
        
    }
    
    void evalNstars(char type, unsigned long long int gSize = 0) { // switch between scaffold, contig, gap while computing N* statistics
        
        switch(type) {
                
            case 's': {
                
                computeNstars(scaffLens, scaffNstars, scaffLstars, &scaffNGstars, &scaffLGstars, gSize);
                break;
                
            }
                
            case 'c': {
                
                computeNstars(contigLens, contigNstars, contigLstars, &contigNGstars, &contigLGstars, gSize);
                break;
                
            }
                
            case 'g': {
                
                computeNstars(gapLens, gapNstars, gapLstars);
                break;
                
            }
                
        }
        
    }
    
    
    void computeNstars(std::vector<unsigned int>& lens, // compute N/L* statistics, vector of all lengths
                       std::vector<unsigned int>& Nstars,      std::vector<unsigned int>& Lstars, // required arguments are passed by reference
                       std::vector<unsigned int>* NGstars = 0, std::vector<unsigned int>* LGstars = 0, unsigned long long int gSize = 0) { // optional arguments are passed by pointer
        
        sort(lens.begin(), lens.end(), std::greater<unsigned int>()); // sort lengths Z-A
        
        unsigned long long int sum = 0, totLen = 0;
        
        for(std::vector<unsigned int>::iterator it = lens.begin(); it != lens.end(); ++it) // find total length
            totLen += *it;
        
        short int N = 1, NG = 1;
        
        for(unsigned int i = 0; i < lens.size(); i++) { // for each length
            
            sum += lens[i]; // increase sum
            
            while (sum >= ((double) totLen / 10 * N) && N<= 10) { // conditionally add length.at or pos to each N/L* bin
                
                Nstars[N-1] = lens[i];
                Lstars[N-1] = i + 1;
                
                N = N + 1;
                
            }
            
            while (gSize > 0 && (sum >= ((double) gSize / 10 * NG)) && NG<= 10) { // if not computing gap statistics repeat also for NG/LG* statistics
                
                (*NGstars)[NG-1] = lens[i];
                (*LGstars)[NG-1] = i + 1;
                
                NG = NG + 1;
                
            }
            
        }
        
    }
    
    void evalAuN(char type, unsigned long long int gSize = 0) { // switch between scaffold, contig, gap while computing N* statistics
        
        switch(type) {
                
            case 's': {
                
                computeAuN(scaffLens, scaffAuN, &scaffAuNG, gSize);
                break;
                
            }
                
            case 'c': {
                
                computeAuN(contigLens, contigAuN, &contigAuNG, gSize);
                break;
                
            }
                
            case 'g': {
                
                computeAuN(gapLens, gapAuN);
                break;
                
            }
                
        }
        
    }
    
    void computeAuN(std::vector<unsigned int>& lens, double& auN, double* auNG = 0, unsigned long long int gSize = 0) {// compute N* statistics
        
        unsigned long long int totLen = 0;
        
        for(std::vector<unsigned int>::iterator it = lens.begin(); it != lens.end(); ++it) // find total length
            totLen += *it;
        
        for(unsigned int i = 0; i < lens.size(); i++) { // for each length
            
            auN += (double) lens[i] * lens[i] / totLen; // the area under the curve is the length (height) * fraction of the total length
            
            if (gSize > 0) {
                
                *auNG += (double) lens[i] * lens[i] / gSize;
                
            }
            
        }
        
    }
    
    unsigned int getScaffN() {
        
        return scaffN;
        
    }
    
    std::vector <unsigned int> getScaffNstars() {
        
        return scaffNstars;
        
    }
    
    std::vector <unsigned int> getScaffNGstars() {
        
        return scaffNGstars;
        
    }
    
    std::vector <unsigned int> getScaffLstars() {
        
        return scaffLstars;
        
    }
    
    std::vector <unsigned int> getScaffLGstars() {
        
        return scaffLGstars;
        
    }
    
    std::vector <unsigned int> getContigNstars() {
        
        return contigNstars;
        
    }
    
    std::vector <unsigned int> getContigNGstars() {
        
        return contigNGstars;
        
    }
    
    std::vector <unsigned int> getContigLstars() {
        
        return contigLstars;
        
    }
    
    std::vector <unsigned int> getContigLGstars() {
        
        return contigLGstars;
        
    }
    
    std::vector <unsigned int> getGapNstars() {
        
        return gapNstars;
        
    }
    
    std::vector <unsigned int> getGapLstars() {
        
        return gapLstars;
        
    }
    
    unsigned int getScaffN50() {
        
        return scaffNstars[4];
        
    }
    
    unsigned int getScaffNG50() {
        
        return scaffNGstars[4];
        
    }
    
    unsigned int getScaffL50() {
        
        return scaffLstars[4];
        
    }
    
    unsigned int getScaffLG50() {
        
        return scaffLGstars[4];
        
    }
    
    double getScaffauN() {
        
        return scaffAuN;
        
    }
    
    double getScaffauNG() {
        
        return scaffAuNG;
        
    }
    
    double getContigauN() {
        
        return contigAuN;
        
    }
    
    double getContigauNG() {
        
        return contigAuNG;
        
    }
    
    double getGapauN() {
        
        return gapAuN;
        
    }
    
    unsigned int getContigNs() {
        
        return inSegments.size();
        
    }
    
    unsigned int getContigN50() {
        
        return contigNstars[4]; // middle value
        
    }
    
    unsigned int getContigNG50() {
        
        return contigNGstars[4];
        
    }
    
    unsigned int getContigL50() {
        
        return contigLstars[4];
        
    }
    
    unsigned int getContigLG50() {
        
        return contigLGstars[4];
        
    }
    
    unsigned int getGapN50() {
        
        return gapNstars[4];
        
    }
    
    unsigned int getGapL50() {
        
        return gapLstars[4];
        
    }
    
    unsigned int getLargestScaffold() {
        
        return scaffLens.size() == 0 ? 0 : scaffLens[0]; // sorted during N/L* computation
        
    }
    
    unsigned int getLargestContig() {
        
        return contigLens.size() == 0 ? 0 : contigLens[0]; // sorted during N/L* computation
        
    }
    
    unsigned int getLargestGap() {
        
        return gapLens.size() > 0 ? gapLens[0] : 0; // sorted during N/L* computation
        
    }
    
    double computeAvgScaffLen() {
        
        return (double) totScaffLen/scaffN;
        
    }
    
    long unsigned int getTotContigLen () {
        
        long unsigned int totContigLen = 0;
        
        for (std::vector<unsigned int>::iterator contigLen = contigLens.begin(); contigLen != contigLens.end(); contigLen++) {
            
            totContigLen += *contigLen;
            
        }
        
        return totContigLen;
        
    }

    double computeAvgContigLen() {
        
        return (double) getTotContigLen()/contigLens.size();
        
    }
    
    double computeAvgSegmentLen() {
        
        return (double) totSegmentLen/inSegments.size();
        
    }
    
    unsigned int getTotGapLen() {
        
        long unsigned int totGapLen = 0;
        
        for (std::vector<unsigned int>::iterator gapLen = gapLens.begin(); gapLen != gapLens.end(); gapLen++) {
            
            totGapLen += *gapLen;
            
        }
        
        return totGapLen;
        
    }
    
    double computeAverageGapLen() {
        
        return totGapLen == 0 ? 0 : (double) getTotGapLen()/gapLens.size();
        
    }
    
    unsigned long int getTotA() {
        
        return totA;
    }
    
    unsigned long int getTotC() {
        
        return totC;
    }
    
    unsigned long int getTotG() {
        
        return totG;
    }
    
    unsigned long int getTotT() {
        
        return totT;
    }
    
    unsigned long int getTotLowerCount() {
        
        return totLowerCount;
    }
    
    double computeGCcontent() {
        
        double GCcontent;
        
        if (inSegments.size()>0) {
            
            GCcontent = (double) (totC + totG) / (totA + totC + totG + totT) * 100;
            
        }else{
            
            GCcontent = 0;
            
        }
        
        return GCcontent;
    }
    
    //gfa methods
    bool addGap(InGap inGap) {
        
        recordGapLen(inGap.dist);
        
        verbose("Recorded length of gaps in sequence");
        
        changeTotGapLen(inGap.dist);
        
        verbose("Increased total gap length");
        
        inGaps.push_back(inGap);

        verbose("Gap added to gap vector");
        
        return true;
        
    }
    
    bool addPath(InPath path) {
        
        inPaths.push_back(path);

        verbose("Path added to path vector");
        
        pathN++;
        
        verbose("Increased path counter");
        
        return true;
        
    }
    
    std::vector<InGap> getGaps() { // return gfa gaps vector
        
        return inGaps;
        
    }

    std::vector<InEdge> getEdges() { // return gfa edge vector
        
        return inEdges;
        
    }
    
    bool appendEdge(InEdge edge) {
        
        inEdges.push_back(edge);

        verbose("Edge added to edge vector");
        
        return true;
        
    }
    
    //sorting methods
    
    void sortPathsByNameAscending(){
        
        sort(inPaths.begin(), inPaths.end(), [](InPath& one, InPath& two){return one.getHeader() < two.getHeader();});
        
    }
    
    void sortPathsByNameDescending(){
        
        sort(inPaths.begin(), inPaths.end(), [](InPath& one, InPath& two){return one.getHeader() > two.getHeader();});
        
    }
    
    void sortPathsByList(std::vector<std::string> headerList){
        
        int index1 = 0, index2 = 0;
        
        auto comp = [&](InPath& one, InPath& two)-> bool { // lambda function for custom sorting
        
        auto it = find(headerList.begin(), headerList.end(), one.getHeader());
      
        if (it != headerList.end()) { // if element one was found

            index1 = it - headerList.begin(); // calculating the index

        }else {

            std::cout<<"Error: sequence missing from sort list (" << one.getHeader() << ")\n";
            exit(1);

        }
            
        it = find(headerList.begin(), headerList.end(), two.getHeader());
      
        if (it != headerList.end()) { // if element two was found

            index2 = it - headerList.begin(); // calculating the index

        }else {

            std::cout<<"Error: sequence missing from sort list ("<<two.getHeader()<<")\n";
            exit(1);

        }

            return index1<index2;
            
        };
        
        sort(inPaths.begin(), inPaths.end(), comp);
        
    }
    
    void sortPathsBySize(bool largest){
        
        auto comp = [&](InPath& one, InPath& two)-> bool { // lambda function for custom sorting
        
            std::vector<PathTuple> pathComponents;
            
            unsigned int uId = 0, sIdx = 0, gIdx = 0, size1 = 0, size2 = 0;
                
            pathComponents = one.getComponents();
            
            for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                
                uId = std::get<1>(*component);
                
                if (std::get<0>(*component) == 'S') {
                
                    auto sId = find_if(inSegments.begin(), inSegments.end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                    
                    if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index
                    
                    size1 += inSegments[sIdx].getInSequence().size();
                    
                }else{
                    
                    auto gId = find_if(inGaps.begin(), inGaps.end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                    
                    if (gId != inGaps.end()) {gIdx = std::distance(inGaps.begin(), gId);} // gives us the segment index
                    
                    size1 += inGaps[gIdx].getDist();
                    
                }
                
            }
                
            pathComponents = two.getComponents();
            
            for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                
                uId = std::get<1>(*component);
                
                if (std::get<0>(*component) == 'S') {
                
                    auto sId = find_if(inSegments.begin(), inSegments.end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                    
                    if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index
                    
                    size2 += inSegments[sIdx].getInSequence().size();
                    
                }else{
                    
                    auto gId = find_if(inGaps.begin(), inGaps.end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                    
                    if (gId != inGaps.end()) {gIdx = std::distance(inGaps.begin(), gId);} // gives us the segment index
                    
                    size2 += inGaps[gIdx].getDist();
                    
                }
            
            
            }
                
            if(largest) {
            
                return size1<size2;
                
            }else{
                
                return size1>size2;
                
            }
                
        };
            
        sort(inPaths.begin(), inPaths.end(), comp);
        
    }
    
    //gfa methods
    void insertHash1(std::string segHeader, unsigned int i) {

        headersToIds.insert({segHeader, i});

    }

    void insertHash2(unsigned int i, std::string segHeader) {

        idsToHeaders.insert({i, segHeader});

    }
    
    void setuId(unsigned int uid) {

        uId = uid;

    }
    
    unsigned int getuId() {

        return uId;

    }
    
    phmap::flat_hash_map<std::string, unsigned int> getHash1() {

        return headersToIds;

    }

    phmap::flat_hash_map<unsigned int, std::string> getHash2() {

        return idsToHeaders;

    }
    
    void buildGraph(std::vector<InGap> const& gaps) { // graph constructor
        
        verbose("Started graph construction");
        
        adjListFW.clear();
        adjListBW.clear();
        
        adjListFW.resize(uId); // resize the adjaciency list to hold all nodes
        adjListBW.resize(uId); // resize the adjaciency list to hold all nodes
        
        for (auto &gap: gaps) // add edges to the graph
        {
            
            verbose("Adding forward gap " + std::to_string(gap.uId) + ": " + idsToHeaders[gap.sId1] + "(" + std::to_string(gap.sId1) + ") " + gap.sId1Or + " " + idsToHeaders[gap.sId2] + "(" + std::to_string(gap.sId2) + ") " + gap.sId2Or + " " + std::to_string(gap.dist));
            
            adjListFW.at(gap.sId1).push_back(std::make_tuple(gap.sId1Or, gap.sId2, gap.sId2Or, gap.dist, gap.uId)); // insert at gap start gap destination, orientations and weight (gap size)

            verbose("Adding reverse gap " + std::to_string(gap.uId) + ": " + idsToHeaders[gap.sId2] + "(" + std::to_string(gap.sId2) + ") " + edge.sId2Or + " " + idsToHeaders[gap.sId1] + "(" + std::to_string(gap.sId1) + ") " + edge.sId2Or + " " + std::to_string(gap.dist));
            
            adjListBW.at(gap.sId2).push_back(std::make_tuple(gap.sId2Or, gap.sId1, gap.sId1Or, gap.dist, gap.uId)); // undirected graph
            
        }
        
        verbose("Graph built");
        
        visited.clear();
        
    }

    void buildEdgeGraph(std::vector<InEdge> const& edges) { // graph constructor
        
        verbose("Started edge graph construction");
        
        adjEdgeListFW.clear();
        
        adjEdgeListFW.resize(uId); // resize the adjaciency list to hold all nodes
        
        for (auto &edge: edges) // add edges to the graph
        {
            
            verbose("Adding forward edge: " + idsToHeaders[edge.sId1] + "(" + std::to_string(edge.sId1) + ") " + edge.sId1Or + " " + idsToHeaders[edge.sId2] + "(" + std::to_string(edge.sId2) + ") " + edge.sId2Or);
            
            adjEdgeListFW.at(edge.sId1).push_back(std::make_tuple(edge.sId1Or, edge.sId2, edge.sId2Or)); // insert at edge start gap destination and orientations
            
        }
        
        verbose("Graph built");
        
        visited.clear();
        
//        assignIds(); // this is not used at present, and seems outdated (segfault on some templates)
        
    }

   void dfsEdges(unsigned int v, unsigned int* componentLength) { // Depth First Search to explore graph connectivity

       visited[v] = true; // mark the current node as visited

       unsigned int sIdx = 0;

       auto sId = find_if(inSegments.begin(), inSegments.end(), [v](InSegment& obj) {return obj.getuId() == v;}); // given a node Uid, find it
        
       if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index

       if (adjEdgeListFW.at(v).size() > 1) { // if the vertex has more than one edge

            *componentLength += inSegments[sIdx].getSegmentLen(); 

            char sign = std::get<0>(adjEdgeListFW.at(v).at(0));
            unsigned int i = 0;

            for(EdgeTuple edge : adjEdgeListFW.at(v)) {
            
                i++;

                if(std::get<0>(edge) != sign){

                    verbose("node: " + idsToHeaders[v] + " --> case a: internal node, multiple edges");

                    break;

                }else if (i == adjEdgeListFW.at(v).size()) {

                    verbose("node: " + idsToHeaders[v] + " --> case b: single dead end, multiple edges");

                    deadEnds += 1;
                }

            sign = std::get<0>(edge);

            }

        }else if (adjEdgeListFW.at(v).size() == 1){ // this is a single dead end

            deadEnds += 1;
            *componentLength += inSegments[sIdx].getSegmentLen(); 

            verbose("node: " + idsToHeaders[v] + " --> case c: single dead end, single edge");

        }else if(adjEdgeListFW.at(v).size() == 0){ // disconnected component (double dead end)

            deadEnds += 2;

            disconnectedComponents++;
            lengthDisconnectedComponents += inSegments[sIdx].getSegmentLen(); 

            verbose("node: " + idsToHeaders[v] + " --> case d: disconnected component");

        }

        for (EdgeTuple i: adjEdgeListFW[v]) { // recur for all forward vertices adjacent to this vertex

           if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {

               dfsEdges(std::get<1>(i), componentLength); // recurse

           }
        }
    }
    
    void dfsScaffolds(unsigned int v, unsigned int* scaffSize, unsigned int* A, unsigned int* C, unsigned int* G, unsigned int* T, unsigned int* lowerCount) // Depth First Search to explore graph connectivity
    {
        
        visited[v] = true; // mark the current node as visited
        unsigned int idx = 0, a = 0, c = 0, g = 0, t = 0;
        
        bool seqRevCom = false, segRevCom = false;
        
        auto it = find_if(inSegments.begin(), inSegments.end(), [&v](InSegment& obj) {return obj.getuId() == v;}); // given a vertex id, search it in the segment vector
        
        if (it != inSegments.end()) {idx = std::distance(inSegments.begin(), it);} // if found, get its index
        
        if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && !backward) { // if the vertex has exactly one forward and one backward connection and they do not connect to the same vertex (internal node)
            
            verbose("node: " + idsToHeaders[v] + " --> case a: internal node, forward direction");
            
            seqRevCom = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? false : true; // check if sequence should be in forward orientation, if not reverse-complement
            
            if (seqRevCom) {
                
                unsigned int tmpA = *A, tmpC = *C;
                
                *A = *T;
                *C = *G;
                *G = tmpC;
                *T = tmpA;
                
            }
            
            segRevCom = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? false : true; // check if vertex should be in forward orientation, if not reverse-complement
            
            backward = false;
            
        }else if (adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 1){ // this is the final vertex without gaps
            
            verbose("node: " + idsToHeaders[v] + " --> case b: end node, forward direction, no final gap");
            
            segRevCom = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? false : true;
            
            backward = true; // reached the end
            
        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 2){ // this is the final vertex with terminal gap
            
            verbose("node: " + idsToHeaders[v] + " --> case c: end node, forward direction, final gap");
            
            if (std::get<1>(adjListBW.at(v).at(0)) != v) { // make sure you are not using the terminal edge to ascertain direction in case it was edited by sak
            
                segRevCom = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? false : true;
             
                *scaffSize += std::get<3>(adjListBW.at(v).at(1));
                
            }else{
            
                segRevCom = (std::get<0>(adjListBW.at(v).at(1)) == '+') ? false : true;
            
                *scaffSize += std::get<3>(adjListBW.at(v).at(0));
                
            }
            
            backward = true; // reached the end
            
        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && backward){ // this is an intermediate vertex, only walking back
            
            verbose("node: " + idsToHeaders[v] + " --> case d: intermediate node, backward direction");
            
            segRevCom = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? false : true;
            
            backward = true;
            
        }else if(adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 0){ // disconnected component
            
            verbose("node: " + idsToHeaders[v] + " --> case e: disconnected component");
            
        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 0){ // this is the first vertex without gaps
            
            verbose("node: " + idsToHeaders[v] + " --> case f: start node, no gaps");
            
            segRevCom = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? false : true;
            
            backward = false;
            
        }else if (adjListFW.at(v).size() == 2 && adjListBW.at(v).size() == 1){ // this is the first vertex with a terminal gap
            
            verbose("node: " + idsToHeaders[v] + " --> case g: start node, start gap");
            
            segRevCom = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? false : true;
            
            *scaffSize += std::get<3>(adjListFW.at(v).at(0));
            
            backward = false;
            
        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) { // if the vertex has exactly one forward and one backward connection and they connect to the same vertex (disconnected component with gap)
            
            verbose("node: " + idsToHeaders[v] + " --> case h: disconnected component with gap");
            
            segRevCom = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? false : true;
            
            *scaffSize += std::get<3>(adjListFW.at(v).at(0));
            
            backward = false;
            
        }
        
        *scaffSize += inSegments[idx].getInSequence().size();
        
        if (!segRevCom) {
            
            a = inSegments[idx].getA();
            c = inSegments[idx].getC();
            g = inSegments[idx].getG();
            t = inSegments[idx].getT();
            
        }else {
            
            a = inSegments[idx].getT();
            c = inSegments[idx].getG();
            g = inSegments[idx].getC();
            t = inSegments[idx].getA();
            
        }
        
        *A += a;
        *C += c;
        *G += g;
        *T += t;
        
        *lowerCount += inSegments[idx].getLowerCount();
        
        for (Tuple i: adjListFW[v]) { // recur for all forward vertices adjacent to this vertex
            
            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {
                
                *scaffSize += std::get<3>(i);
                
                dfsScaffolds(std::get<1>(i), scaffSize, A, C, G, T, lowerCount); // recurse
                
            }
        }
        
        for (Tuple i: adjListBW[v]) { // recur for all backward vertices adjacent to this vertex
            
            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {
                
                *scaffSize += std::get<3>(i);
                
                dfsScaffolds(std::get<1>(i), scaffSize, A, C, G, T, lowerCount); // recurse
                
            }
        }
        
    }
    
    std::vector<std::vector<Tuple>> getAdjListFW() {
        
        return adjListFW;
        
    }
    
    std::vector<std::vector<Tuple>> getAdjListBW() {
        
        return adjListBW;
        
    }
    
    bool getVisited(unsigned int uId) {
        
        return visited[uId];
        
    }
    
    bool getDeleted(unsigned int uId) {
        
        return deleted[uId];
        
    }
    
    bool updateScaffoldStats() {
        
        scaffLens.clear();
        contigLens.clear();
        gapLens.clear();
        
        totScaffLen = 0, scaffN = 0, totA = 0, totC = 0, totG = 0, totT = 0, totLowerCount = 0;
        
        for (InPath& inPath : inPaths) { // loop through all paths
            
            walkPath(inPath.getpUId());
            
            totScaffLen += inPath.getLen();
         
            verbose("Increased total scaffold length");
            
            scaffN++;
            
            verbose("Increased total scaffold N");
            
            recordScaffLen(inPath.getLen());
            
            verbose("Recorded length of sequence: " + std::to_string(inPath.getLen()));
            
            totA += inPath.getA();
            totC += inPath.getC();
            totG += inPath.getG();
            totT += inPath.getT();
            
            verbose("Increased total ACGT counts");
            
            totLowerCount += inPath.getLowerCount();

            verbose("Increased total count of lower bases");
            
        }
        
        verbose("Updated scaffold statistics");
        
        return true;
        
    }
    
    bool updateGapLens() {
        
        totGapLen = 0;
        gapLens.clear();
        
        for (unsigned int i = 0; i != inGaps.size(); ++i) { // loop through all edges
        
            recordGapLen(inGaps[i].getDist());
            
            changeTotGapLen(inGaps[i].getDist());
            
        }
        
        verbose("Updated list of gap lengths");
        
        return true;
        
    }
    
    bool removeTerminalGaps() { // if two contigs are provided, remove all edges connecting them, if only one contig is provided remove all edges where it appears
        
        std::vector<InPath>::iterator pathIt = inPaths.begin(); // first, remove the gaps from the paths they occur in
        std::vector<PathTuple> pathComponents;
        
        unsigned int uId = 0, gIdx = 0;
        
        while (pathIt != end(inPaths)) {
            
            pathComponents = pathIt->getComponents();
            
            for (std::vector<PathTuple>::iterator componentIt = pathComponents.begin(); componentIt != pathComponents.end();) {
                
                if (std::get<0>(*componentIt) == 'G') {
                    
                    uId = std::get<1>(*componentIt);
                    
                    auto gId = find_if(inGaps.begin(), inGaps.end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a gap Uid, find it
                    
                    if (gId->getsId1() == gId->getsId2()) {
                        
                        gIdx = std::distance(pathComponents.begin(), componentIt); // gives us the gap index
                    
                        pathIt->pathComponents.erase(pathIt->pathComponents.begin()+gIdx); // remove the gap component from the path
                        
                        verbose("Removed gap from paths");
                        
                        inGaps.erase(gId); // remove the element by position, considering elements that were already removed in the loop

                        changeTotGapLen(-gId->getDist()); // update length of gaps
                        
                        verbose("Removed gap from gap vector");
                        
                    }
                    
                }
                
                componentIt++;
                
            }
            
            pathIt++;
            
        }
        
        gapLens.clear();
        
        for (unsigned int i = 0; i != inGaps.size(); ++i) { // update statistics to reflect the removal
        
            recordGapLen(inGaps[i].getDist());
            
        }
        
        verbose("Recorded length of gaps in sequence");
        
        return true;
        
    }

    unsigned int getDeadEnds() {

        return deadEnds;

    }

    unsigned int getDisconnectedComponents() {

        return disconnectedComponents;

    }
    
    unsigned int getLengthDisconnectedComponents() {

        return lengthDisconnectedComponents;

    }

    // instruction methods
    
    std::vector<InGap> getGap(std::string* contig1, std::string* contig2 = NULL) { // if two contigs are provided, returns all edges connecting them, if only one contig is provided returns all edges where it appears
 
        std::vector<InGap> gaps;
        
        unsigned int sUId1 = headersToIds[*contig1];
        
        if (contig2 != NULL) {
        
            unsigned int sUId2 = headersToIds[*contig2];
        
            auto gId = find_if(inGaps.begin(), inGaps.end(), [sUId1, sUId2](InGap& obj) {return ( // given two vertex ids, search the gap that connects them
                
                (obj.getsId1() == sUId1 && obj.getsId2() == sUId2) || // fw orientation
                (obj.getsId1() == sUId2 && obj.getsId2() == sUId1)    // rv orientation
                                                                                                 
            );});
            
            gaps.push_back(*gId);
            
        }else{
            
            auto it = inGaps.begin();
            
            while (it != end(inGaps)) {

                auto gId = find_if(it, inGaps.end(), [sUId1](InGap& obj) {return obj.getsId1() == sUId1 || obj.getsId2() == sUId1;}); // check whether an edge containing the node was found
                
                if (gId != inGaps.end()) {
                
                    gaps.push_back(*gId);
                    
                }
                
                it++;
                
            }
            
        }
        
        return gaps;
        
    }
    
    std::vector<unsigned int> removeGaps(std::string* contig1, std::string* contig2 = NULL) { // if two contigs are provided, remove all edges connecting them, if only one contig is provided remove all edges where it appears
        
        std::vector<unsigned int> guIds;
 
        unsigned int sUId1 = headersToIds[*contig1], gIdx = 0;
        
        if (contig2 != NULL) {
        
            unsigned int sUId2 = headersToIds[*contig2];
        
            auto gId = find_if(inGaps.begin(), inGaps.end(), [sUId1, sUId2](InGap& obj) {return ( // given two vertex ids, search the gap that connects them
                
                (obj.getsId1() == sUId1 && obj.getsId2() == sUId2) || // fw orientation
                (obj.getsId1() == sUId2 && obj.getsId2() == sUId1)    // rv orientation
                                                                                                 
            );});
        
            if (gId != inGaps.end()) {
                
                gIdx = std::distance(inGaps.begin(), gId); // gives us the gap index
                
                guIds.push_back((*gId).getuId());
            
                inGaps.erase(inGaps.begin()+gIdx); // remove the element by position, considering elements that were already removed in the loop
                
                changeTotGapLen(-gId->getDist()); // update length of gaps
                
            }
            
        }else{
            
            auto it = inGaps.begin();
            
            while (it != end(inGaps)) {

                auto gId = find_if(it, inGaps.end(), [sUId1](InGap& obj) {return obj.getsId1() == sUId1 || obj.getsId2() == sUId1;}); // check whether an edge containing the node was found

                if (gId != inGaps.end()) {
                    
                    gIdx = std::distance(inGaps.begin(), gId); // gives us the gap index
                    
                    guIds.push_back((*gId).getuId());
            
                    inGaps.erase(inGaps.begin()+gIdx); // remove the element by position, considering elements that were already removed in the loop
                    
                    changeTotGapLen(-gId->getDist()); // update length of gaps
                    
                }
                
                it = gId;
                
            }
            
        }
        
        return guIds;
        
    }
    
    bool removeSegment(std::string* contig1) { // if two contigs are provided, remove all edges connecting them, if only one contig is provided remove all edges where it appears
        
        if (contig1 != NULL) {
            
            unsigned int sIdx = 0, sUId = headersToIds[*contig1];
        
            auto sId = find_if(inSegments.begin(), inSegments.end(), [sUId](InSegment& obj) {return obj.getuId() == sUId;}); // given a node uId, find it
        
            if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index
        
            deleted[sIdx] = true;
            
            changeTotSegmentLen(-sId->getSegmentLen());
            
        }else{
            
            verbose("Cannot detect node: " + *contig1);
            
        }
        
        return true;
        
    }
    
    void removePath(unsigned int pUId, bool silent = false) {
        
        auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
        
        if (pathIt != inPaths.end()) {
            
            inPaths.erase(pathIt);
            
        }else if (!silent){
            
            fprintf(stderr, "Warning: the path you are attempting to remove does not exist (pUId: %i). Skipping.\n", pUId);
            
        }
    
    }
    
    void removeGap(unsigned int gUId, bool silent = false) {
        
        auto gapIt = find_if(inGaps.begin(), inGaps.end(), [gUId](InGap& obj) {return obj.getuId() == gUId;}); // given a path pUId, find it
        
        if (gapIt != inGaps.end()) {
            
            inGaps.erase(gapIt);
            
        }else if (!silent){
            
            fprintf(stderr, "Warning: the gap you are attempting to remove does not exist (pUId: %i). Skipping.\n", gUId);
            
        }
    
    }
    
    void removePathsFromSegment(unsigned int uId) {
        
        int i = 0;
        
        for (InPath inPath : inPaths) {
            
            std::vector<PathTuple> pathComponents = inPath.getComponents();
            
            auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [uId](PathTuple& obj) {return std::get<1>(obj) == uId;}); // given a node uId, find if present in the given path
            
            if (pathIt != pathComponents.end()) {
            
                removePath(i);
                
            }
            
            ++i;
            
        }
        
    }
    
    void removeSegmentInPath(unsigned int uId, InGap gap) {
        
        int i = 0;
        
        InPath newPath;
        
        std::vector<PathTuple> newComponents;
        
        for (InPath inPath : inPaths) {
            
            std::vector<PathTuple> pathComponents = inPath.getComponents();
            
            auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [uId](PathTuple& obj) {return std::get<1>(obj) == uId;}); // given a node uId, find if present in the given path
            
            if (pathIt != pathComponents.end()) {
                
                newPath.setHeader(inPath.getHeader()); // set path header
                
                newPath.append({std::begin(pathComponents), pathIt-1}); // subset path excluding segment to be removed
                
                newPath.add('G', gap.getuId(), '0'); // introduce gap
                
                newPath.append({pathIt+2, std::end(pathComponents)}); // add remaining components
                
                addPath(newPath);
                
                removePath(i);
                
                break;
                    
            }
            
            ++i;
            
        }
        
    }
    
    void joinPaths(std::string pHeader, unsigned int pUId1, unsigned int pUId2, std::string gHeader, unsigned int gUId, char pId1Or, char pId2Or, unsigned int dist, unsigned int start1, unsigned int end1, unsigned int start2, unsigned int end2) {
        
        InPath path;
        
        phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find (pHeader); // get the headers to uIds table to look for the header
        
        if (got == headersToIds.end()) { // this is the first time we see this path
            
            verbose("Path not found in keys. Creating new path (" + pHeader + ", pUId: " + std::to_string(uId) + ")");
            
            insertHash1(pHeader, uId); // header to hash table
            insertHash2(uId, pHeader); // header to hash table
            
            path.newPath(uId, pHeader);
            
            uId++;
            
        }else{
            
            path.setHeader(pHeader);
            pUId1 = got->second;
            path.setpUId(pUId1);
            
            verbose("Path already exists in keys (" + pHeader + ", pUId: " + std::to_string(pUId1) + "). Joining.");
            
        }
        
        PathTuple component1, component2;
            
        auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId1](InPath& obj) {return obj.getpUId() == pUId1;}); // given a path pUId, find it
        
        if (pathIt != inPaths.end()) {
            
            verbose("Path found in path set (pIUd: " + std::to_string(pUId1) + "). Adding components to new path.");
            
            std::vector<PathTuple> pathComponents = pathIt->getComponents();
            
            trimPathByRef(pathComponents, start1, end1);
            
            if (pId1Or == '-') {revComPathComponents(pathComponents);}
            
            path.append({std::begin(pathComponents), std::end(pathComponents)});
            
            component1 = *std::prev(std::end(pathComponents));
            
            if (start1 == 0) {
            
                removePath(pUId1); // remove path1
                
            }
            
        }else{
            
            fprintf(stderr, "Warning: could not locate in path set (pIUd: %u)\n", pUId1);
            
        }
        
        if (gUId == 0) {
            
            insertHash1(gHeader, uId); // header to hash table
            insertHash2(uId, gHeader); // id to hash table
            
            gUId = uId;
            
            uId++;
        
        }
        
        verbose("Adding gap to new path (" + gHeader + ")");
        
        path.add('G', gUId, '0');
            
        pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId2](InPath& obj) {return obj.getpUId() == pUId2;}); // given a path pUId, find it
        
        if (pathIt != inPaths.end()) {
            
            verbose("Path found in path set (pIUd: " + std::to_string(pUId2) + "). Adding components to new path.");
            
            std::vector<PathTuple> pathComponents = pathIt->getComponents();
            
            trimPathByRef(pathComponents, start2, end2);
            
            if (pId2Or == '-') {revComPathComponents(pathComponents);}
        
            path.append({std::begin(pathComponents), std::end(pathComponents)});
            
            component2 = *std::begin(pathComponents);
            
            if (start2 == 0) {
            
                removePath(pUId2); // remove path2
                
            }
            
        }else{
            
            fprintf(stderr, "Warning: could not locate in path set (pIUd: %u)\n", pUId2);
            
        }
        
        InGap gap;
        
        gap.newGap(gUId, std::get<1>(component1), std::get<1>(component2), '+', '+', dist, gHeader); // define the new gap
        
        addGap(gap); // introduce the new gap
        
        addPath(path);
        
    }
    
    InPath joinPathsByComponent(std::string seqHeader, unsigned int uId1, unsigned int uId2, unsigned int uId3) {
        
        int i = 0;
        
        InPath newPath;
        newPath.setHeader(seqHeader);
        
        for (InPath inPath : inPaths) { // add first path
            
            std::vector<PathTuple> pathComponents = inPath.getComponents();
            
            auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [uId1](PathTuple& obj) {return std::get<1>(obj) == uId1;}); // given a node uId, find if present in the given path
            
            if (pathIt != pathComponents.end()) {
            
                newPath.append({std::begin(pathComponents), std::end(pathComponents)});
            
                break;
                
            }
            
            ++i;
            
        }
        
        newPath.add('G', uId2, '0');
        
        for (InPath inPath : inPaths) { // add second path
            
            std::vector<PathTuple> pathComponents = inPath.getComponents();
            
            auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [uId3](PathTuple& obj) {return std::get<1>(obj) == uId3;}); // given a node uId, find if present in the given path
            
            if (pathIt != pathComponents.end()) {
            
                newPath.append({std::begin(pathComponents), std::end(pathComponents)});
            
                break;
                
            }
            
            ++i;
            
        }
        
        return newPath;
        
    }
    
    void splitPath(unsigned int guId, std::string header1, std::string header2) {
        
        int i = 0;
        
        InPath newPath1;
        newPath1.setHeader(header1);
        std::vector<PathTuple> newComponents1;
        
        InPath newPath2;
        newPath2.setHeader(header2);
        std::vector<PathTuple> newComponents2;
        
        for (InPath inPath : inPaths) { // add first path
            
            std::vector<PathTuple> pathComponents = inPath.getComponents();
            
            auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [guId](PathTuple& obj) {return std::get<1>(obj) == guId;}); // given a node uId, find if present in the given path
            
            if (pathIt != pathComponents.end()) {
            
                newComponents1.insert(std::end(newComponents1), std::begin(pathComponents), pathIt);
                
                newPath1.setComponents(newComponents1);
                
                if (newComponents1.size() > 0) {
                
                    addPath(newPath1);
                    
                }
                
                newComponents2.insert(std::end(newComponents2), pathIt+1, std::end(pathComponents));
                
                newPath2.setComponents(newComponents2);
                
                if (newComponents2.size() > 0) {
                
                    addPath(newPath2);
                    
                }
                
                removePath(i);
            
                break;
                
            }
            
            ++i;
            
        }
        
    }
    
    void clearPaths() {
        
        inPaths.clear();
        
    }
    
    void clearGaps() {
        
        inGaps.clear();
        
    }
    
    void renamePath(unsigned int pUId, std::string pHeader, unsigned int* newpUId) {
        
        auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
        
        pathIt->setHeader(pHeader);
        
        insertHash1(pHeader, pUId); // header to hash table
        insertHash2(pUId, pHeader); // uId to hash table
        
        if (newpUId != NULL) {
            
            pathIt->setpUId(*newpUId);
            
        }
        
    }
    
    void revComPath(unsigned int pUId) {
        
        auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
        
        pathIt->revCom();
        
        verbose("Path reverse-complemented (" + pathIt->getHeader() + ").");
        
    }
    
    void trimPathByUId(unsigned int pUId, unsigned int start, unsigned int end) {
        
        auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
        
        std::vector<PathTuple>* pathComponents = pathIt->getComponentsByRef();
        
        trimPath(pathComponents, start, end);
        
    }

    void trimPathByRef(std::vector<PathTuple>& pathComponents, unsigned int start, unsigned int end) {
        
        trimPath(&pathComponents, start, end);
        
    }

    void trimPath(std::vector<PathTuple>* pathComponents, unsigned int start, unsigned int end) {
        
        if(start == 0 || end == 0) {

            verbose("Nothing to trim. Skipping.");
            return;

        }
        
        verbose("Trimming path (start: " + std::to_string(start) + ", end: " + std::to_string(end) + ")");
    
        unsigned int cUId = 0, traversedSize = 0, actualSize = 0, compSize = 0, newCompSize = 0, compOriginalSize = 0;
        
        for (std::vector<PathTuple>::iterator component = pathComponents->begin(); component != pathComponents->end(); component++) {
            
            cUId = std::get<1>(*component);
            
            verbose("New path size before iteration: " + std::to_string(actualSize));
            
            verbose("Checking original coordinates of component (uId: " + std::to_string(std::get<1>(*component)) + ", start: " + std::to_string(std::get<3>(*component)) + ", end: " + std::to_string(std::get<4>(*component)) + ")");
            
            if (std::get<0>(*component) == 'S') {
            
                auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                compOriginalSize = inSegment->getSegmentLength();
                
                compSize = inSegment->getSegmentLength(std::get<3>(*component), std::get<4>(*component));
                
                if (std::get<4>(*component) > 0) {
                    
                    verbose("Component was already trimmed (size: " + std::to_string(compSize) + ")");
                    
                }else{
                    
                    verbose("Component was not already trimmed (size: " + std::to_string(compSize) + ")");
                    
                }
                
            }else{
                
                auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                    
                compSize = inGap->getDist(std::get<3>(*component), std::get<4>(*component));
                
                compOriginalSize = inGap->getDist();
                
            }
                
            if (traversedSize + compSize < start) {
                
                verbose("Start coordinate exceeds component, removing it");
                
                pathComponents->erase(component);
                
                component--;
                
                traversedSize += compSize;
                
                continue;
                
            }

            if (traversedSize + compSize >= start && traversedSize < start) {
                
                if (std::get<2>(*component) == '+') { // we only change the end coordinate if the component wasn't already flipped, otherwise we edit the start
                
                std::get<3>(*component) = (std::get<3>(*component) == 0 ? 0 : std::get<3>(*component) - 1) + start - traversedSize; // edit also to account for already trimmed component
                    
                }else{
                    
                    std::get<4>(*component) = end - traversedSize;
                    
                }
                
                verbose("Start coordinate of the component needs to be edited as result of subsetting (new start: " + std::to_string(std::get<3>(*component)) + ")");
                
                newCompSize = std::get<4>(*component) - std::get<3>(*component) + 1;
                
            }
            
            if (std::get<3>(*component) > 0 && std::get<4>(*component) == 0) { // account for editing of the end coordinate but not the start coordinate
                
                std::get<4>(*component) = compSize;
                
            }
            
            if (traversedSize + compSize >= end) {
                
                pathComponents->erase(component + 1, pathComponents->end());
                
                verbose("Erased extra components");
                
                if (std::get<2>(*component) == '+') { // we only change the end coordinate if the component wasn't already flipped, otherwise we edit the start
                
                    std::get<4>(*component) = end - traversedSize; // edit also to account for already trimmed component // traversedSize < start ? 0 : std::get<3>(*component)
                    
                }else{
                    
                    std::get<3>(*component) = (std::get<3>(*component) == 0 ? 0 : std::get<3>(*component) - 1) + start - traversedSize; // edit also to account for already trimmed component
                    
                }
                
                if (std::get<3>(*component) > std::get<4>(*component)) {std::get<4>(*component) = std::get<3>(*component) + std::get<4>(*component);}
                    
                verbose("End coordinate of the component needs to be edited as result of subsetting (new end: " + std::to_string(std::get<4>(*component)) + ")");
                
                if (std::get<4>(*component) > 0 && std::get<3>(*component) == 0) { // account for editing of the end coordinate but not the start coordinate
                    
                    std::get<3>(*component) = 1;
                    
                }
                
                newCompSize = std::get<4>(*component) - std::get<3>(*component) + 1;
                
            }
            
            if(std::get<3>(*component) == 1 && std::get<4>(*component) == compOriginalSize) { // if the result of trimming restores the original size of the component, no need to adjust the internal coordinates
                
                std::get<3>(*component) = 0;
                std::get<4>(*component) = 0;
                
            }
            
            actualSize += newCompSize;
            verbose("Path size after iteration: " + std::to_string(actualSize));
            
            traversedSize += compSize;
            verbose("Traversed path: " + std::to_string(traversedSize));
            
        }
            
        verbose("Final path size: " + std::to_string(actualSize));
    
    }
    
    unsigned int pathLen(unsigned int pUId) {
        
        auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
        
        return pathIt->getLen();
        
    }
    
    void walkPath(unsigned int pUId) {
        
        unsigned int cUId = 0, gapLen = 0;

        auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
        
        std::vector<PathTuple> pathComponents = pathIt->getComponents();
        
        for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                
            cUId = std::get<1>(*component);
        
            if (std::get<0>(*component) == 'S') {
                
                auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                contigLens.push_back(inSegment->getSegmentLength(std::get<3>(*component), std::get<4>(*component)));
                
                pathIt->increaseLen(inSegment->getSegmentLength(std::get<3>(*component), std::get<4>(*component)));
                
                pathIt->increaseLowerCount(inSegment->getLowerCount(std::get<4>(*component) - std::get<3>(*component)));
                
                if (std::get<3>(*component) == 0 || std::get<4>(*component) == 0) {
                
                    pathIt->increaseA(std::get<2>(*component) == '+' ? inSegment->getA() : inSegment->getT());
                    
                    pathIt->increaseC(std::get<2>(*component) == '+' ? inSegment->getC() : inSegment->getG());
                    
                    pathIt->increaseG(std::get<2>(*component) == '+' ? inSegment->getG() : inSegment->getC());
                    
                    pathIt->increaseT(std::get<2>(*component) == '+' ? inSegment->getT() : inSegment->getA());
                    
                }else{
                    
                    std::string sequence = inSegment->getInSequence(std::get<3>(*component), std::get<4>(*component));
                    
                    if (std::get<2>(*component) == '+') {
                    
                        for (char base : sequence) {
                            
                            switch (base) {
                                case 'A':
                                case 'a':{
                                    
                                    pathIt->increaseA(1);
                                    break;
                                    
                                }
                                case 'C':
                                case 'c':{
                                    
                                    pathIt->increaseC(1);
                                    break;
                                    
                                }
                                case 'G':
                                case 'g': {
                                    
                                    pathIt->increaseG(1);
                                    break;
                                    
                                }
                                case 'T':
                                case 't': {
                                    
                                    pathIt->increaseT(1);
                                    break;
                                    
                                }
                        
                            }
                            
                        }

                    }else{
                        
                        for (char base : sequence) {
                        
                            switch (base) {
                                case 'A':
                                case 'a':{
                                    
                                    pathIt->increaseT(1);
                                    break;
                                    
                                }
                                case 'C':
                                case 'c':{
                                    
                                    pathIt->increaseG(1);
                                    break;
                                    
                                }
                                case 'G':
                                case 'g': {
                                    
                                    pathIt->increaseC(1);
                                    break;
                                    
                                }
                                case 'T':
                                case 't': {
                                    
                                    pathIt->increaseA(1);
                                    break;
                                    
                                }
                        
                            }
                            
                        }
                        
                    }
                    
                }
                
            }else{
                
                auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                gapLen += inGap->getDist(std::get<4>(*component) - std::get<3>(*component));
                
                if (!(std::get<0>(*(component + 1)) == 'G')) {
                
                    gapLens.push_back(gapLen);
                    
                    gapLen = 0;
                    
                }
                
                pathIt->increaseLen(inGap->getDist(std::get<4>(*component) - std::get<3>(*component)));
                
            }
            
        }
        
    }
    
    void discoverPaths() {
        
        buildGraph(inGaps);
        
        for (InSegment inSegment : inSegments) {
            
            if (!visited[inSegment.getuId()]) {
                
                InPath path;
                
                path.newPath(uId, inSegment.getSeqHeader() + "_path");
                
                insertHash1(inSegment.getSeqHeader() + "_path", uId); // header to hash table
                insertHash2(uId, inSegment.getSeqHeader() + "_path"); // uId to hash table
                
                uId++;
                
                dfsPath(inSegment.getuId(), path);
                
                addPath(path);
                
            }
            
        }
        
    }
    
    void dfsPath(unsigned int v, InPath& newPath) // Depth First Search to build a new path given a vertex
    {

        visited[v] = true; // mark the current node as visited

        if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && !backward) { // if the vertex has exactly one forward and one backward connection and they do not connect to the same vertex (internal node)

            verbose("node: " + idsToHeaders[v] + " --> case a: internal node, forward direction");
            
            newPath.add('S', v, '+');

            backward = false;

        }else if (adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 1){ // this is the final vertex without gaps

            verbose("node: " + idsToHeaders[v] + " --> case b: end node, forward direction, no final gap");
            
            newPath.add('S', v, '-');

            backward = true; // reached the end

        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 2){ // this is the final vertex with terminal gap

            verbose("node: " + idsToHeaders[v] + " --> case c: end node, forward direction, final gap");
            
            newPath.add('S', v, '-');

            if (std::get<1>(adjListBW.at(v).at(0)) != v) { // make sure you are not using the terminal edge to ascertain direction in case it was edited by sak

            }else{

            }

            backward = true; // reached the end

        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && backward){ // this is an intermediate vertex, only walking back

            verbose("node: " + idsToHeaders[v] + " --> case d: intermediate node, backward direction");
            
            newPath.add('S', v, '-');

            backward = true;

        }else if(adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 0){ // disconnected component
            
            verbose("node: " + idsToHeaders[v] + " --> case e: disconnected component");
            
            newPath.add('S', v, '+');

        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 0){ // this is the first vertex without gaps
            
            verbose("node: " + idsToHeaders[v] + " --> case f: start node, no gaps");
            
            newPath.add('S', v, '+');

            visited.clear();

            visited[v] = true; // we have just visited the start node

            backward = false;

        }else if (adjListFW.at(v).size() == 2 && adjListBW.at(v).size() == 1){ // this is the first vertex with a terminal gap

            verbose("node: " + idsToHeaders[v] + " --> case g: start node, start gap");

            newPath.add('S', v, '+');

            visited.clear();

            visited[v] = true; // we have just visited the start node

            backward = false;

        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) { // if the vertex has exactly one forward and one backward connection and they connect to the same vertex (disconnected component with gap)

            verbose("node: " + idsToHeaders[v] + " --> case h: disconnected component with gap");
            
            newPath.add('S', v, '+');

            backward = false;

        }

        for (Tuple i: adjListFW[v]) { // recur for all forward vertices adjacent to this vertex

            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {

                dfsPath(std::get<1>(i), newPath); // recurse

            }
        }

        for (Tuple i: adjListBW[v]) { // recur for all backward vertices adjacent to this vertex

            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {

                dfsPath(std::get<1>(i), newPath); // recurse

            }
        }

    }
    
    void findBubbles () {
        
        unsigned int sUId = 0, sUId1 = 0, sUId2 = 0;
        char sId1Or, sId2Or;
        
        visited.clear();
        
        buildEdgeGraph(inEdges);
        
        for (InSegment segment : inSegments) {
            
            sUId = segment.getuId();
            
            verbose("Evaluating for node: " + idsToHeaders[sUId] + " (uId: " + std::to_string(sUId) + ")");
            
            if (!visited[sUId] && !deleted[sUId]) { // check if the node was already visited
                
                verbose("The node was not yet visited or deleted");
                
                if (adjEdgeListFW.at(sUId).size() == 2 && std::get<0>(adjEdgeListFW.at(sUId).at(0)) != std::get<0>(adjEdgeListFW.at(sUId).at(1))) { // if it has exactly two edges with different orientation it could be a het region of the bubble
                    
                    verbose("Exactly two edges with different orientation found. Could be a het region of the bubble");
                    
                    // then check the the adjacient nodes
                    sUId1 = std::get<1>(adjEdgeListFW.at(sUId).at(0));
                    sId1Or = std::get<2>(adjEdgeListFW.at(sUId).at(0));
                    
                    sUId2 = std::get<1>(adjEdgeListFW.at(sUId).at(1));
                    sId2Or = std::get<2>(adjEdgeListFW.at(sUId).at(1));
                    
                    if (adjEdgeListFW.at(sUId1).size() >= 2 && adjEdgeListFW.at(sUId2).size() >= 2) { // both nodes need at least two edges to be a bubble
                        
                        verbose("Both neighbour nodes have at least two edges");
                        
                        for (EdgeTuple edge1 : adjEdgeListFW.at(sUId1)) {
                            
                            if (std::get<2>(edge1) == sId1Or && std::get<1>(edge1) != sUId) { // we are checking edges on the side of the potential bubble for node1, avoiding the original node
                                
                                verbose("Evaluating node: " + idsToHeaders[std::get<1>(edge1)] + " (uId: " + std::to_string(std::get<1>(edge1)) + ")");
                                
                                if (std::get<1>(edge1) == sUId2) { // this is a potential insertion
                                    
                                    bubbles.push_back(std::make_tuple(sUId, sUId1, sUId2, NULL));
                                    
                                    verbose("Candidate insertion found");
                                    
                                }else{
                                    
                                    for (EdgeTuple edge2 : adjEdgeListFW.at(sUId2)) {
                                    
                                        if (std::get<2>(edge2) == sId2Or && std::get<1>(edge1) == std::get<1>(edge2)) { // we are checking edges on the side of the potential bubble for node2 and that it connects to the same node as node1
                                            
                                            verbose("Evaluating node: " + idsToHeaders[std::get<1>(edge2)] + " (uId: " + std::to_string(std::get<1>(edge2)) + ")");
                                            
                                            bubbles.push_back(std::make_tuple(sUId, sUId1, sUId2, std::get<1>(edge2)));
                                            
                                            verbose("Candidate bubble found");
                                            
                                        }
                                        
                                    }
                                    
                                }
                                   
                            }
                            
                        }
                        
                    }
                    
                }

            }
            
            visited[sUId] = true;
            
        }
        
    }
    
    std::vector<Bubble>* getBubbles () {
        
        return &bubbles;
        
    }
    
    // end of gfa methods
    
};

#endif /* GFASTATS_GFA_H */
