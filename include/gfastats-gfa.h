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
    
    std::string getInSequence() {
        return inSequence;
    }
    
    std::string getInSequenceQuality() {
        
        return inSequenceQuality;
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
    
    unsigned int getLowerCount() {
        
        return lowerCount;
    }
    
    unsigned int getSegmentLength() {
        
        return inSequence.size();
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
    
    unsigned int getDist() {
        
        return dist;
        
    }
    
};
class InEdge {
    private:
//    unsigned long long int lineN; // useful if we wish to sort as is the original input
    std::string cigar;
    char sId1Or, sId2Or;
    unsigned int eUId, eId, sId1, sId2;
    
    friend class SAK;
    friend class InSequences;
    
public:
    void newEdge(unsigned int eUid, unsigned int sid1, unsigned int sid2, const char& sid1or, const char& sid2or, std::string c = "") {
        
        eUId = eUid;
        sId1 = sid1;
        sId2 = sid2;
        sId1Or = sid1or;
        sId2Or = sid2or;
        cigar = c;
        
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

public:
    
    void setHeader(std::string pheader) {
        
        pHeader = pheader;
    
    }
    
    void setComment(std::string c) {
        pComment = c;
    }
    
    void addToPath(char type, unsigned int UId, char sign) {
        
        pathComponents.push_back(std::make_tuple(type, UId, sign));
        
    }
    
    void clearPath() {
        
        pathComponents.clear();
        
    }
    
    std::vector<PathTuple> getComponents() {
        
        return pathComponents;
        
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
    std::unordered_map<std::string, unsigned int> headersToIds;
    std::unordered_map<unsigned int, std::string> idsToHeaders;
    std::unordered_map<int, bool> visited, deleted;
    bool backward = false, first = false;
    
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
    
    unsigned long long int totSegmentLen = 0;
    
    unsigned int
    scaffN = 0,
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
        
        inSegment.setLowerCount(lowerCount);
        
        verbose("Increased count of lower bases");
        
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
        
        path.addToPath('G', *uId, '+');
        
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
        
        path.addToPath('S', *uId, '+');
        
        *A = 0, *C = 0, *G = 0, *T = 0;
        
        (*iId)++; // number of gaps in the current scaffold
        (*uId)++; // unique numeric identifier
        
    }
    
    void traverseInSequence(std::string* seqHeader, std::string* seqComment, std::string* sequence, std::string* sequenceQuality = NULL) { // traverse the sequence to split at gaps and measure sequence properties
        
        unsigned int pos = 0, // current position in sequence
        A = 0, C = 0, G = 0, T = 0,
        lowerCount = 0,
        dist = 0, // gap size
        iId = 1, // scaffold feature internal identifier
        sStart = 0, sEnd = 0; // segment start and end
        char sign = '+';
        bool wasN = false;
        
        path.clearPath();
        path.setHeader(*seqHeader);

        if (seqComment != NULL) {
            
            path.setComment(*seqComment);
            
        }
        
        unsigned int seqLen = sequence->length()-1;
        
        for (char &base : *sequence) {
            
            if (islower(base)) {
                
                lowerCount++;
                
            }
            
            switch (base) {
                    
                case 'N':
                case 'n':
                case 'X':
                case 'x': {
                    
                    dist++;
                    
                    if (!wasN && pos>0) { // gap start and gap not at the start of the sequence
                            
                        sEnd = pos - 1;
                        pushbackSegment(seqHeader, seqComment, sequence, &iId, &uId, &A, &C, &G, &T, &lowerCount, sStart, sEnd, sequenceQuality);
                        
                    }
                    
                    if(pos == seqLen) { // end of scaffold, terminal gap
                        
                        sign = '-';
                        
                        pushbackGap(seqHeader, &iId, &uId, pos, &dist, sign, seqLen, -1);
                        
                    }
                    
                    wasN = true;
                    
                    break;
                }
                default: {
                    
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
                    
                    if (wasN) { // internal gap end
                        
                        sStart = pos;
                        pushbackGap(seqHeader, &iId, &uId, pos, &dist, sign, seqLen, 1);
                        
                    }
                    
                    if (pos == seqLen) {
                        
                        sEnd = pos;
                        pushbackSegment(seqHeader, seqComment, sequence, &iId, &uId, &A, &C, &G, &T, &lowerCount, sStart, sEnd, sequenceQuality);
                        
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
        
        std::unordered_map<std::string, unsigned int>::const_iterator got = headersToIds.find (*seqHeader); // get the headers to uIds table a look for the header
        
        if (got == headersToIds.end()) { // this is the first time we see this segment
            
            insertHash1(*seqHeader, uId); // header to hash table
            insertHash2(uId, *seqHeader); // header to hash table
            
            sUId = uId;
            
            uId++;
            
        }else{
            
            sUId = got->second;
            
        }
                
        addSegment(sUId, 0, *seqHeader, seqComment, sequence, &A, &C, &G, &T, &lowerCount, sequenceQuality);
            
        A = 0, C = 0, G = 0, T = 0;
        
    }
    
    void appendSequence(std::string* seqHeader, std::string* seqComment, std::string* sequence, std::string* sequenceQuality = NULL) { // method to append a new sequence from a fasta
        
        verbose("Header, comment, sequence and (optionally) quality read");
        
        if(verbose_flag) {std::cout<<"\n";};
        
        traverseInSequence(seqHeader, seqComment, sequence, sequenceQuality);
        
        verbose("Sequence traversed");
        
        scaffN++;
        
        verbose("Increased total scaffold N");
        
        if(verbose_flag) {std::cout<<"\n";};
        
    }
    
    void appendSegment(std::string* seqHeader, std::string* seqComment, std::string* sequence, std::string* sequenceQuality = NULL) { // method to append a new segment from a gfa
        
        verbose("Header, comment, sequence and (optionally) quality read");
        
        if(verbose_flag) {std::cout<<"\n";};
        
        traverseInSegment(seqHeader, seqComment, sequence, sequenceQuality);
        
        verbose("Segment traversed");
        
        if(verbose_flag) {std::cout<<"\n";};
        
    }
    
    unsigned int getSegmentN() {
        
        return inSegments.size();
    }
    
    InSegment getInSegment(unsigned int idx) {
        
        InSegment inSegment = inSegments[idx];
        return inSegment;
        
    }
    
    std::vector<InSegment> getInSegments() {
        
        return inSegments;
        
    }
    
    std::vector<InGap> getInGaps() {
        
        return inGaps;
        
    }
    
    std::vector<InPath> getInPaths() {
        
        return inPaths;
        
    }
    
    unsigned long long int getTotScaffLen() {
        
        return totSegmentLen+totGapLen;
        
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
    
    unsigned int getTotGapLen() {
        
        return totGapLen;
        
    }
    
    unsigned int getGapN() {
        
        return inGaps.size();
        
    }

    unsigned int getEdgeN() {
        
        return inEdges.size();
        
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
        
        return scaffLens[0]; // sorted during N/L* computation
        
    }
    
    unsigned int getLargestContig() {
        
        return contigLens[0]; // sorted during N/L* computation
        
    }
    
    unsigned int getLargestGap() {
        
        return gapLens.size() > 0 ? gapLens[0] : 0; // sorted during N/L* computation
        
    }
    
    double computeAverageScaffLen() {
        
        return (double) (totSegmentLen+totGapLen)/scaffN;
        
    }
    
    double computeAverageSegmentLen() {
        
        return (double) totSegmentLen/contigLens.size();
        
    }
    
    double computeAverageGapLen() {
        
        return (double) totGapLen/gapLens.size();
        
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
        
        verbose("Added nodes to hash table");
        
        recordGapLen(inGap.dist);
        
        verbose("Recorded length of gaps in sequence");
        
        changeTotGapLen(inGap.dist);
        
        verbose("Increased total gap length");
        
        inGaps.push_back(inGap);

        verbose("Gap added to gap vector");
        
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
    
    std::unordered_map<std::string, unsigned int> getHash1() {

        return headersToIds;

    }

    std::unordered_map<unsigned int, std::string> getHash2() {

        return idsToHeaders;

    }
    
    void buildGraph(std::vector<InGap> const& edges) { // graph constructor
        
        verbose("Started graph construction");
        
        adjListFW.clear();
        adjListBW.clear();
        
        adjListFW.resize(inSegments.size()+inGaps.size()+inEdges.size()); // resize the adjaciency list to hold all nodes
        adjListBW.resize(inSegments.size()+inGaps.size()+inEdges.size()); // resize the adjaciency list to hold all nodes
        
        for (auto &edge: edges) // add edges to the graph
        {
            
            verbose("Adding forward edge: " + idsToHeaders[edge.sId1] + "(" + std::to_string(edge.sId1) + ") " + edge.sId1Or + " " + idsToHeaders[edge.sId2] + "(" + std::to_string(edge.sId2) + ") " + edge.sId2Or + " " + std::to_string(edge.dist));
            
            adjListFW.at(edge.sId1).push_back(std::make_tuple(edge.sId1Or, edge.sId2, edge.sId2Or, edge.dist)); // insert at gap start gap destination, orientations and weight (gap size)

            verbose("Adding reverse edge: " + idsToHeaders[edge.sId2] + "(" + std::to_string(edge.sId2) + ") " + edge.sId2Or + " " + idsToHeaders[edge.sId1] + "(" + std::to_string(edge.sId1) + ") " + edge.sId2Or + " " + std::to_string(edge.dist));
            
            adjListBW.at(edge.sId2).push_back(std::make_tuple(edge.sId2Or, edge.sId1, edge.sId1Or, edge.dist)); // undirected graph
            
        }
        
        verbose("Graph built");
        
        visited.clear();
        
//        assignIds(); // this is not used at present, and seems outdated (segfault on some templates)
        
    }
    
//    void assignIds() { // (re)assignment, whenever the graph is built or modified
//
//        unsigned int i = 0, uId = 0;
//
//        for (InSegment inSegment : inSegments) {// loop through all nodes
//
//            if (!visited[i] && !deleted[i]) { // if the node was not visited yet
//
//                dfsPos(i, &uId); // start from the first node and try to assign sId = 0
//
//                uId = 0;
//
//            }
//
//            i++;
//
//        }
//
//        backward = false;
//
//        visited.clear();
//
//        verbose("Assigned node IDs");
//
//    }

    void buildEdgeGraph(std::vector<InEdge> const& edges) { // graph constructor
        
        verbose("Started edge graph construction");
        
        adjEdgeListFW.clear();
        
        adjEdgeListFW.resize(inSegments.size()+inGaps.size()+inEdges.size()); // resize the adjaciency list to hold all nodes
        
        for (auto &edge: edges) // add edges to the graph
        {
            
            verbose("Adding forward edge: " + idsToHeaders[edge.sId1] + "(" + std::to_string(edge.sId1) + ") " + edge.sId1Or + " " + idsToHeaders[edge.sId2] + "(" + std::to_string(edge.sId2) + ") " + edge.sId2Or);
            
            adjEdgeListFW.at(edge.sId1).push_back(std::make_tuple(edge.sId1Or, edge.sId2, edge.sId2Or)); // insert at edge start gap destination and orientations
            
        }
        
        verbose("Graph built");
        
        visited.clear();
        
//        assignIds(); // this is not used at present, and seems outdated (segfault on some templates)
        
    }

   void dfsEdges(int v) { // Depth First Search to find graph connectivity

       visited[v] = true; // mark the current node as visited

       if (adjEdgeListFW.at(v).size() > 1) { // if the vertex has more than one edge

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

           verbose("node: " + idsToHeaders[v] + " --> case c: single dead end, single edge");

       }else if(adjEdgeListFW.at(v).size() == 0){ // disconnected component (double dead end)

            deadEnds += 2;

            disconnectedComponents++;

            verbose("node: " + idsToHeaders[v] + " --> case d: disconnected component");

       }

       for (EdgeTuple i: adjEdgeListFW[v]) { // recur for all forward vertices adjacent to this vertex

           if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {

               dfsEdges(std::get<1>(i)); // recurse

           }
       }

   }
    
//    void dfsPos(int v, unsigned int* uId) { // Depth First Search to assign sequential sIds to each feature in the graph. It finds the first node then walks to the end node assigning progressive uIds
//
//        visited[v] = true; // mark the current node as visited
//
//        (*uId)++; // increase the segment sId counter
//
//        if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && !backward) { // if the vertex has exactly one forward and one backward connection and they do not connect to the same vertex (internal node)
//
//            verbose("node: " + idsToHeaders[v] + " --> case a: internal node, forward direction");
//
//            backward = false;
//
//        }else if (adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 1){ // this is the final vertex without gaps
//
//            verbose("node: " + idsToHeaders[v] + " --> case b: end node, forward direction, no final gap");
//
//            backward = true; // reached the end
//
//        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 2){ // this is the final vertex with terminal gap
//
//            verbose("node: " + idsToHeaders[v] + " --> case c: end node, forward direction, final gap");
//
//            backward = true; // reached the end
//
//        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && backward){ // this is an intermediate vertex, only walking back
//
//            verbose("node: " + idsToHeaders[v] + " --> case d: intermediate node, backward direction");
//
//            backward = true;
//
//        }else if(adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 0){ // disconnected component
//
//            verbose("node: " + idsToHeaders[v] + " --> case e: disconnected component");
//
//        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 0){ // this is the first vertex without gaps
//
//            verbose("node: " + idsToHeaders[v] + " --> case f: start node, no gaps");
//
//            visited.clear();
//
//            *uId = 1; // we have just found the start node, count features from here
//
//            visited[v] = true; // we have just visited the start node
//
//            backward = false;
//
//        }else if (adjListFW.at(v).size() == 2 && adjListBW.at(v).size() == 1){ // this is the first vertex with a terminal gap
//
//            verbose("node: " + idsToHeaders[v] + " --> case g: start node, start gap");
//
//            visited.clear();
//
//            *uId = 1; // we have just found the start node, count nodes from here
//
//            visited[v] = true; // we have just visited the start node
//
//            backward = false;
//
//        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) { // if the vertex has exactly one forward and one backward connection and they connect to the same vertex (disconnected component with gap)
//
//            verbose("node: " + idsToHeaders[v] + " --> case h: disconnected component with gap");
//
//            backward = false;
//
//        }
//
//        inSegments[v].setuId(*uId); // set the uId for the node
//
//        verbose("Assigned node internal id: " + std::to_string(*uId));
//
//        for (Tuple i: adjListFW[v]) { // recur for all forward vertices adjacent to this vertex
//
//            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {
//
//                inGaps[v].setuId(*uId); // set the uId for the edge
//
//                (*uId)++; // increase the edge Id counter
//
//                dfsPos(std::get<1>(i), uId); // recurse
//
//            }
//        }
//
//        for (Tuple i: adjListBW[v]) { // recur for all backward vertices adjacent to this vertex
//
//            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {
//
//                dfsPos(std::get<1>(i), uId); // recurse
//
//            }
//        }
//
//    }
    
    
    // LEGACY, we use paths now (still potentially useful if path is absent
    
//    void dfsSeq(unsigned int v, std::string &inSequence, std::string* inSequenceQuality = NULL) // Depth First Search to build fast* sequence
//    {
//
//        visited[v] = true; // mark the current node as visited
//        std::string inSequenceNext;
//        std::string inSequenceQualityNext;
//        unsigned int idx = 0;
//
//        auto it = find_if(inSegments.begin(), inSegments.end(), [&v](InSegment& obj) {return obj.getuId() == v;}); // given a vertex id, search it in the segment vector
//
//        if (it != inSegments.end()) {idx = std::distance(inSegments.begin(), it);} // if found, get its index
//
//        if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && !backward) { // if the vertex has exactly one forward and one backward connection and they do not connect to the same vertex (internal node)
//
//            verbose("node: " + idsToHeaders[v] + " --> case a: internal node, forward direction");
//
//            inSequence = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? inSequence : revCom(inSequence); // check if sequence should be in forward orientation, if not reverse-complement
//
//            inSequenceNext = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence()); // check if vertex should be in forward orientation, if not reverse-complement
//
//            inSequence += inSequenceNext;
//
//            if (!(inSequenceQuality == NULL)) {
//
//                *inSequenceQuality = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? *inSequenceQuality : rev(*inSequenceQuality); // check if vertex should be in forward orientation, if not reverse-complement
//
//                inSequenceQualityNext = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? inSegments[idx].getInSequenceQuality() : rev(inSegments[idx].getInSequenceQuality()); // check if vertex should be in forward orientation, if not reverse-complement
//
//                *inSequenceQuality += inSequenceQualityNext;
//
//            }
//
//            backward = false;
//
//        }else if (adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 1){ // this is the final vertex without gaps
//
//            verbose("node: " + idsToHeaders[v] + " --> case b: end node, forward direction, no final gap");
//
//            inSequenceNext = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence());
//
//            inSequence += inSequenceNext;
//
//            if (!(inSequenceQuality == NULL)) {
//
//                inSequenceQualityNext = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? inSegments[idx].getInSequenceQuality() : rev(inSegments[idx].getInSequenceQuality());
//
//                *inSequenceQuality += inSequenceQualityNext;
//
//            }
//
//            backward = true; // reached the end
//
//        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 2){ // this is the final vertex with terminal gap
//
//            verbose("node: " + idsToHeaders[v] + " --> case c: end node, forward direction, final gap");
//
//            if (std::get<1>(adjListBW.at(v).at(0)) != v) { // make sure you are not using the terminal edge to ascertain direction in case it was edited by sak
//
//                inSequenceNext = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence());
//
//                inSequence += inSequenceNext;
//
//                inSequence += std::string(std::get<3>(adjListBW.at(v).at(1)), 'N'); // add gap
//
//            }else{
//
//                inSequenceNext = (std::get<0>(adjListBW.at(v).at(1)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence());
//
//                inSequence += inSequenceNext;
//
//                inSequence += std::string(std::get<3>(adjListBW.at(v).at(0)), 'N'); // add gap
//
//            }
//
//            if (!(inSequenceQuality == NULL)) {
//
//                if (std::get<1>(adjListBW.at(v).at(0)) != v) { // make sure you are not using the terminal edge to ascertain direction in case it was edited by sak
//
//                    inSequenceQualityNext = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? inSegments[idx].getInSequenceQuality() : rev(inSegments[idx].getInSequenceQuality());
//
//                    *inSequenceQuality += inSequenceQualityNext;
//
//                    *inSequenceQuality += std::string(std::get<3>(adjListBW.at(v).at(1)), '!'); // add missing quality
//
//                }else{
//
//                    inSequenceQualityNext = (std::get<0>(adjListBW.at(v).at(1)) == '+') ? inSegments[idx].getInSequenceQuality() : rev(inSegments[idx].getInSequenceQuality());
//
//                    *inSequenceQuality += inSequenceQualityNext;
//
//                    *inSequenceQuality += std::string(std::get<3>(adjListBW.at(v).at(0)), '!'); // add missing quality
//
//                }
//
//            }
//
//            backward = true; // reached the end
//
//        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && backward){ // this is an intermediate vertex, only walking back
//
//            verbose("node: " + idsToHeaders[v] + " --> case d: intermediate node, backward direction");
//
//            inSequenceNext = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence());
//
//            inSequence.insert(0, inSequenceNext);
//
//            if (!(inSequenceQuality == NULL)) {
//
//                inSequenceQualityNext = (std::get<0>(adjListBW.at(v).at(0)) == '+') ? inSegments[idx].getInSequenceQuality() : rev(inSegments[idx].getInSequenceQuality());
//
//                (*inSequenceQuality).insert(0, inSequenceQualityNext);
//
//            }
//
//            backward = true;
//
//        }else if(adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 0){ // disconnected component
//
//            verbose("node: " + idsToHeaders[v] + " --> case e: disconnected component");
//
//            inSequence += inSegments[idx].getInSequence();
//
//            if (!(inSequenceQuality == NULL)) {
//
//                *inSequenceQuality += inSegments[idx].getInSequenceQuality();
//
//            }
//
//        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 0){ // this is the first vertex without gaps
//
//            verbose("node: " + idsToHeaders[v] + " --> case f: start node, no gaps");
//
//            inSequenceNext = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence());
//
//            inSequence.insert(0, inSequenceNext);
//
//            if (!(inSequenceQuality == NULL)) {
//
//                inSequenceQualityNext = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? inSegments[idx].getInSequenceQuality() : rev(inSegments[idx].getInSequenceQuality());
//
//                (*inSequenceQuality).insert(0, inSequenceQualityNext);
//
//            }
//
//            backward = false;
//
//        }else if (adjListFW.at(v).size() == 2 && adjListBW.at(v).size() == 1){ // this is the first vertex with a terminal gap
//
//            verbose("node: " + idsToHeaders[v] + " --> case g: start node, start gap");
//
//            inSequenceNext = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence());
//
//            inSequence.insert(0, inSequenceNext);
//
//            inSequence.insert(0, std::string(std::get<3>(adjListBW.at(v).at(0)), 'N')); // add gap
//
//            if (!(inSequenceQuality == NULL)) {
//
//                (*inSequenceQuality).insert(0, std::string(std::get<3>(adjListBW.at(v).at(0)), '!')); // add missing quality
//
//                inSequenceQualityNext = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? inSegments[idx].getInSequenceQuality() : rev(inSegments[idx].getInSequenceQuality());
//
//                (*inSequenceQuality).insert(0, inSequenceQualityNext);
//
//            }
//
//            backward = false;
//
//        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) { // if the vertex has exactly one forward and one backward connection and they connect to the same vertex (disconnected component with gap)
//
//            verbose("node: " + idsToHeaders[v] + " --> case h: disconnected component with gap");
//
//            inSequenceNext = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence());
//
//            inSequence += inSequenceNext;
//
//            std::get<2>(adjListFW.at(v).at(0)) == '-' ? inSequence += std::string(std::get<3>(adjListFW.at(v).at(0)), 'N') : inSequence.insert(0,std::string(std::get<3>(adjListBW.at(v).at(0)), 'N')); // add gap: if - in second vertex is terminal else is start gap
//
//            if (!(inSequenceQuality == NULL)) {
//
//                inSequenceQualityNext = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? inSegments[idx].getInSequenceQuality() : rev(inSegments[idx].getInSequenceQuality());
//
//                *inSequenceQuality += inSequenceQualityNext;
//
//                std::get<2>(adjListFW.at(v).at(0)) == '-' ? *inSequenceQuality += std::string(std::get<3>(adjListFW.at(v).at(0)), '!') : (*inSequenceQuality).insert(0,std::string(std::get<3>(adjListBW.at(v).at(0)), '!')); // add missing sequence: if - in second vertex is terminal else is start gap
//
//            }
//
//            backward = false;
//
//        }
//
//        for (Tuple i: adjListFW[v]) { // recur for all forward vertices adjacent to this vertex
//
//            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {
//
//                inSequence += std::string(std::get<3>(i), 'N'); // add gaps
//
//                if (!(inSequenceQuality == NULL)) {
//
//                    *inSequenceQuality += std::string(std::get<3>(i), '!'); // add missing quality
//
//                }
//
//                dfsSeq(std::get<1>(i), inSequence, inSequenceQuality); // recurse
//
//            }
//        }
//
//        for (Tuple i: adjListBW[v]) { // recur for all backward vertices adjacent to this vertex
//
//            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {
//
//                inSequence.insert(0, std::string(std::get<3>(i), 'N')); // add gaps
//
//                if (!(inSequenceQuality == NULL)) {
//
//                    (*inSequenceQuality).insert(0, std::string(std::get<3>(i), '!')); // add missing quality
//
//                }
//
//                dfsSeq(std::get<1>(i), inSequence, inSequenceQuality); // recurse
//
//            }
//        }
//
//    }
    
    void dfsAgp(unsigned int v, std::string &agp, unsigned int &cStart, unsigned int &cEnd) // Depth First Search to generate AGP output
    {
        
        visited[v] = true; // mark the current node as visited
        std::string agpNext, seqHeader;
        unsigned int idx = 0;
        
        auto it = find_if(inSegments.begin(), inSegments.end(), [&v](InSegment& obj) {return obj.getuId() == v;}); // given a vertex id, search its index in the segment vector
        
        if (it != inSegments.end()) {idx = std::distance(inSegments.begin(), it);} // gives us the vertex index
        
        seqHeader = inSegments[idx].getSeqHeader();
        
        if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && !backward) { // if the vertex has exactly one forward and one backward connection and they do not connect to the same vertex (internal node)
            
            verbose("node: " + idsToHeaders[v] + " --> case a: internal node, forward direction");
            
            if (first) {
            
                cEnd = cStart + inSegments[idx].getInSequence().size() - 1;
                
                agpNext = seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+inSegments[idx].getSeqHeader().substr(inSegments[idx].getSeqHeader().length() - 1)+"\tW\t"+inSegments[idx].getSeqHeader()+"\t1\t"+std::to_string(inSegments[idx].getInSequence().size())+"\t"+std::string(1, std::get<0>(adjListBW.at(v).at(0)))+"\n";
                
                agp += agpNext;
                
                cStart = cEnd + 1;
                
            }
            
            backward = false;
            
        }else if (adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 1){ // this is the final vertex without gaps
            
            verbose("node: " + idsToHeaders[v] + " --> case b: end node, forward direction, no final gap");
            
            if (first) {
            
                cEnd = cStart + inSegments[idx].getInSequence().size() - 1;
                
                agpNext = seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+inSegments[idx].getSeqHeader().substr(inSegments[idx].getSeqHeader().length() - 1)+"\tW\t"+inSegments[idx].getSeqHeader()+"\t1\t"+std::to_string(inSegments[idx].getInSequence().size())+"\t"+std::string(1, std::get<0>(adjListBW.at(v).at(0)))+"\n";
                
                agp += agpNext;
                
                cStart = cEnd + 1;
            
            }
                
            backward = true; // reached the end
            
        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 2){ // this is the final vertex with terminal gap
            
            verbose("node: " + idsToHeaders[v] + " --> case c: end node, forward direction, final gap");
            
            if (first) {
            
                cEnd = cStart + inSegments[idx].getInSequence().size() - 1;
                
                agpNext = seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+inSegments[idx].getSeqHeader().substr(inSegments[idx].getSeqHeader().length() - 1)+"\tW\t"+inSegments[idx].getSeqHeader()+"\t1\t"+std::to_string(inSegments[idx].getInSequence().size())+"\t"+std::string(1, std::get<2>(adjListBW.at(v).at(0)))+"\n";
                
                agp += agpNext;
                
                cStart = cEnd + 1;
                
                cEnd = cStart + std::get<3>(adjListFW.at(v).at(0)) - 1;
                
                agp += seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+std::to_string(std::get<3>(adjListFW.at(v).at(0)))+"\tN\t"+std::to_string(std::get<3>(adjListFW.at(v).at(0)))+"\tscaffold\tyes\n"; // add gaps
                
                cStart = cEnd + 1;
                
            }
            
            backward = true; // reached the end
            
        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) && backward){ // this is an intermediate vertex, only walking back
            
            verbose("node: " + idsToHeaders[v] + " --> case d: intermediate node, backward direction, doing nothing");
            
            backward = true;
            
        }else if(adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 0){ // disconnected component
            
            verbose("node: " + idsToHeaders[v] + " --> case e: disconnected component");
            
            agpNext = seqHeader+"\t1\t"+std::to_string(inSegments[idx].getInSequence().size())+"\t1\tW\t"+inSegments[idx].getSeqHeader()+"\t1\t"+std::to_string(inSegments[idx].getInSequence().size())+"\t+\n";
            
            agp += agpNext;
            
        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 0){ // this is the first vertex without gaps
            
            verbose("node: " + idsToHeaders[v] + " --> case f: start node, no gaps");
            
            cStart = 1;
            
            cEnd = cStart + inSegments[idx].getInSequence().size() - 1;
            
            agpNext = seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+inSegments[idx].getSeqHeader().substr(inSegments[idx].getSeqHeader().length() - 1)+"\tW\t"+inSegments[idx].getSeqHeader()+"\t1\t"+std::to_string(inSegments[idx].getInSequence().size())+"\t"+std::get<0>(adjListFW.at(v).at(0))+"\n";
            
            agp.insert(0, agpNext);
            
            cStart = cEnd + 1;
            
            backward = false; // we only walk forward now
            
            first = true; // we have identified the first node
            
            visited.clear(); // once the first vertex has been identified restart the walk
            
            visited[v] = true; // we have just visited the start node
            
        }else if (adjListFW.at(v).size() == 2 && adjListBW.at(v).size() == 1){ // this is the first vertex with a start gap
            
            verbose("node: " + idsToHeaders[v] + " --> case g: start node, start gap");
            
            cStart = 1;
            
            cEnd = cStart + std::get<3>(adjListFW.at(v).at(0)) - 1;
            
            agp += seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+std::to_string(std::get<3>(adjListFW.at(v).at(0)))+"\tN\t"+std::to_string(std::get<3>(adjListFW.at(v).at(0)))+"\tscaffold\tyes\n"; // add gaps
            
            cStart = cEnd + 1;
            
            cEnd = cStart + inSegments[idx].getInSequence().size() - 1;
            
            agpNext = (std::get<0>(adjListFW.at(v).at(0)) == '+') ? inSegments[idx].getInSequence() : revCom(inSegments[idx].getInSequence());
            
            agpNext = seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+inSegments[idx].getSeqHeader().substr(inSegments[idx].getSeqHeader().length() - 1)+"\tW\t"+inSegments[idx].getSeqHeader()+"\t1\t"+std::to_string(inSegments[idx].getInSequence().size())+"\t"+std::get<0>(adjListFW.at(v).at(0))+"\n";
            
            agp += agpNext;
            
            cStart = cEnd + 1;
            
            backward = false; // we only walk forward now
            
            first = true; // we have identified the first node
            
            visited.clear(); // once the first vertex has been identified restart the walk
            
            visited[v] = true; // we have just visited the start node
            
        }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && std::get<1>(adjListFW.at(v).at(0)) == std::get<1>(adjListBW.at(v).at(0))) { // if the vertex has exactly one forward and one backward connection and they connect to the same vertex (disconnected component with gap)
            
            verbose("node: " + idsToHeaders[v] + " --> case h: disconnected component with gap");
            
            std::get<2>(adjListFW.at(v).at(0)) == '-' ? // if negative gap (-)
            cStart = 1 :
            cStart = std::get<3>(adjListFW.at(v).at(0)) + 1;
            
            cEnd = cStart + inSegments[idx].getInSequence().size() - 1;
            
            agpNext = seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+inSegments[idx].getSeqHeader().substr(inSegments[idx].getSeqHeader().length() - 1)+"\tW\t"+inSegments[idx].getSeqHeader()+"\t1\t"+std::to_string(inSegments[idx].getInSequence().size())+"\t+\n";
            
            agp += agpNext;
            
            cStart = cEnd + 1;
            
            cEnd = cStart + std::get<3>(adjListFW.at(v).at(0)) - 1;
            
            std::get<2>(adjListFW.at(v).at(0)) == '-' ? // if negative gap (-)
            agp += seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t1\tN\t"+std::to_string(std::get<3>(adjListFW.at(v).at(0)))+"\tscaffold\tyes\n": // insert at the end
            agp.insert(0,seqHeader+"\t1\t"+std::to_string(std::get<3>(adjListFW.at(v).at(0)))+"\t1\tN\t"+std::to_string(std::get<3>(adjListFW.at(v).at(0)))+"\tscaffold\tyes\n"); // insert at the beginning
            
            backward = false;
            
        }
        
        for (Tuple i: adjListFW[v]) { // recur for all forward vertices adjacent to this vertex
            
            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {
                
                if (first) {
                
                    cEnd = cStart + std::get<3>(i) - 1;
                
                    agp += seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+std::to_string(std::get<3>(i))+"\tN\t"+std::to_string(std::get<3>(i))+"\tscaffold\tyes\n"; // add gaps
                
                    cStart = cEnd + 1;
                    
                }
                
                dfsAgp(std::get<1>(i), agp, cStart, cEnd); // recurse
                
            }
        }
        
        for (Tuple i: adjListBW[v]) { // recur for all backward vertices adjacent to this vertex
            
            if (!visited[std::get<1>(i)] && !deleted[std::get<1>(i)]) {
                
                if (first) {
                
                    cEnd = cEnd + std::get<3>(i) - 1;
                
                    agp += seqHeader+"\t"+std::to_string(cStart)+"\t"+std::to_string(cEnd)+"\t"+std::to_string(std::get<3>(i))+"\tN\t"+std::to_string(std::get<3>(i))+"\tscaffold\tyes\n"; // add gaps
                
                    cStart = cEnd + 1;
                    
                }
                
                dfsAgp(std::get<1>(i), agp, cStart, cEnd); // recurse
                
            }
        }
        
    }
    
    void dfsScaffolds(unsigned int v, unsigned int* scaffSize, unsigned int* A, unsigned int* C, unsigned int* G, unsigned int* T, unsigned int* lowerCount) // Depth First Search to build fast* sequence
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
        
        scaffN = 0;
        
        unsigned int scaffSize = 0, A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
        
        buildGraph(getGaps()); // first build the graph
        
        for (InSegment inSegment : inSegments) { // loop through all nodes
            
            if (!getVisited(inSegment.getuId()) && !getDeleted(inSegment.getuId())) { // check if the node was already visited and not deleted
        
                dfsScaffolds(inSegment.getuId(), &scaffSize, &A, &C, &G, &T, &lowerCount); // then walk the scaffold to update statistics
             
                scaffN++;
                
                verbose("Increased total scaffold N");
                
                recordScaffLen(scaffSize);
                
                verbose("Recorded length of sequence: " + std::to_string(scaffSize));
                
                scaffSize = 0;
                
                totA += A;
                totC += C;
                totG += G;
                totT += T;
                
                verbose("Increased ACGT counts");
                
                A = 0, C = 0, G = 0, T = 0;
                
                totLowerCount += lowerCount;

                verbose("Increased count of lower bases");
                
                lowerCount = 0;
                
            }
            
        }
        
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
        
        std::vector<InGap>::iterator it = inGaps.begin();
        
        while (it != end(inGaps)) {
        
            if (it->getsId1() == it->getsId2()) {
            
                inGaps.erase(it); // remove the element by position, considering elements that were already removed in the loop
                
                changeTotGapLen(-it->getDist()); // update length of gaps
                
                it--; // the iterator is now one gap short
                
            }
    
        it++; // check the new gap
            
        }
        
        
        gapLens.clear();
        
        for (unsigned int i = 0; i != inGaps.size(); ++i) { // loop through all edges
        
            recordGapLen(inGaps[i].getDist());
            
            verbose("Recorded length of gaps in sequence");
            
        }
        
        return true;
        
    }

    unsigned int getDeadEnds() {

        return deadEnds;

    }

    unsigned int getDisconnectedComponents() {

        return disconnectedComponents;

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

                auto gId = find_if(it+1, inGaps.end(), [sUId1](InGap& obj) {return obj.getsId1() == sUId1 || obj.getsId2() == sUId1;}); // check whether an edge containing the node was found
                
                if (gId != inGaps.end()) {
                
                    gaps.push_back(*gId);
                    
                }
                
                it = gId;
                
            }
            
        }
        
        return gaps;
        
    }
    
    bool removeGaps(std::string* contig1, std::string* contig2 = NULL) { // if two contigs are provided, remove all edges connecting them, if only one contig is provided remove all edges where it appears
 
        unsigned int sUId1 = headersToIds[*contig1], gIdx = 0;
        
        if (contig2 != NULL) {
        
            unsigned int sUId2 = headersToIds[*contig2];
        
            auto gId = find_if(inGaps.begin(), inGaps.end(), [sUId1, sUId2](InGap& obj) {return ( // given two vertex ids, search the gap that connects them
                
                (obj.getsId1() == sUId1 && obj.getsId2() == sUId2) || // fw orientation
                (obj.getsId1() == sUId2 && obj.getsId2() == sUId1)    // rv orientation
                                                                                                 
            );});
        
            if (gId != inGaps.end()) {
                
                gIdx = std::distance(inGaps.begin(), gId); // gives us the gap index
            
                inGaps.erase(inGaps.begin()+gIdx); // remove the element by position, considering elements that were already removed in the loop
                
                changeTotGapLen(-gId->getDist()); // update length of gaps
                
            }
            
        }else{
            
            auto it = inGaps.begin();
            
            while (it != end(inGaps)) {

                auto gId = find_if(it, inGaps.end(), [sUId1](InGap& obj) {return obj.getsId1() == sUId1 || obj.getsId2() == sUId1;}); // check whether an edge containing the node was found

                if (gId != inGaps.end()) {
                    
                    gIdx = std::distance(inGaps.begin(), gId); // gives us the gap index
            
                    inGaps.erase(inGaps.begin()+gIdx); // remove the element by position, considering elements that were already removed in the loop
                    
                    changeTotGapLen(-gId->getDist()); // update length of gaps
                    
                }
                
                it = gId;
                
            }
            
        }
        
        return true;
        
    }
    
    bool removeSegment(std::string* contig1) { // if two contigs are provided, remove all edges connecting them, if only one contig is provided remove all edges where it appears
        
        if (contig1 != NULL) {
            
            unsigned int sIdx = 0, sUId = headersToIds[*contig1];
        
            auto sId = find_if(inSegments.begin(), inSegments.end(), [sUId](InSegment& obj) {return obj.getuId() == sUId;}); // given a node Uid, find it
        
            if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index
        
            deleted[sIdx] = true;
            
            changeTotSegmentLen(-sId->getSegmentLen());
            
        }else{
            
            verbose("Cannot detect node: " + *contig1);
            
        }
        
        return true;
        
    }
    
    // end of gfa methods
    
};

#endif /* GFASTATS_GFA_H */
