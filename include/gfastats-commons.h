//
//  gfastats-classes.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef gfastatscommons_h
#define gfastatscommons_h

//classes
class BedCoordinates {
private:
    std::vector<std::string> seqHeaders;
    std::vector<unsigned int> cBegin;
    std::vector<unsigned int> cEnd;
    
public:
    
    void pushCoordinates(std::string h, unsigned int b, unsigned int e) {
        
        seqHeaders.push_back(h);
        cBegin.push_back(b);
        cEnd.push_back(e);
        
    }
    
    bool empty() {
        
        return (seqHeaders.size()==0) ? true : false;
        
    }
    
    std::vector<std::string> getSeqHeaders() {
        
        return seqHeaders;
        
    }
    
    std::string getSeqHeader(unsigned int pos) {
        
        return seqHeaders[pos];
        
    }
    
    unsigned int getcBegin(unsigned int pos) {
        
        return cBegin[pos];
        
    }
    
    unsigned int getcEnd(unsigned int pos) {
        
        return cEnd[pos];
        
    }
    
};

class InSequence {
private:
    std::string seqHeader;
    std::string seqComment;
    std::string inSequence;
    std::string inSequenceQuality;
    std::vector<unsigned int> contigBoundaries;
    std::vector<unsigned int> gapBoundaries;
    unsigned int A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
    
public:
    
    void traverseInSequence(std::string* s) {
        
        unsigned int pos = 0, A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
        bool wasN = false, pushbackGap = false;
        std::vector<unsigned int> caseBoundaries;
        std::vector<unsigned int> gapBoundaries;
        gapBoundaries.reserve(200);
            
        for (char &base : *s) {
            
            if (islower(base)) {
                
                lowerCount++;
                
            }
            
            switch (base) {
                    
                case 'N':
                case 'n':
                case 'X':
                case 'x': {
                    
                    if (!wasN) { // gap start
                        
                        pushbackGap = true;
                        
                    }
                    
                    if(pos == (s->length()-1)) { // end of scaffold
                        
                        if (!wasN){
                            
                            gapBoundaries.push_back(pos);
                            
                        }
                        
                        pushbackGap = true;
                        pos++;
                        
                    }
                    
                    wasN = true;
                    
                    break;
                }
                default: {
                    
                    if (wasN) { // internal gap end
                        
                        pushbackGap = true;
                        
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
                    
                    wasN = false;
                    
                }
                    
            }
            
            if (pushbackGap) {
                
                gapBoundaries.push_back(pos);
                pushbackGap = false;
                
            }
            
            pos++;
            
        }
        
        setSeqGapBoundaries(gapBoundaries);
        setSeqContigBoundaries(gapBoundaries);
        setACGT(A, C, G, T);
        setLowerCount(lowerCount);
        
    }
    
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
    
    void setSeqContigBoundaries(std::vector<unsigned int> &gapBoundaries) {
        
        std::vector<unsigned int> newSeqContigBoundaries;
        
        newSeqContigBoundaries.reserve(gapBoundaries.size() + 2);
        
        if (gapBoundaries.size() > 0) {
            
            newSeqContigBoundaries = gapBoundaries;
            
            if (gapBoundaries[0] != 0) {
                
                newSeqContigBoundaries.insert(newSeqContigBoundaries.begin(), 0);
                
            }else{
                
                newSeqContigBoundaries.erase(newSeqContigBoundaries.begin());
                
            }
            
            if (newSeqContigBoundaries[newSeqContigBoundaries.size()-1] != inSequence.size()) {
                
                newSeqContigBoundaries.insert(newSeqContigBoundaries.end(), inSequence.size());
                
            }else{
                
                newSeqContigBoundaries.pop_back();
                
            }
            
        }else{
            
            newSeqContigBoundaries = {0, (unsigned int) inSequence.size()};
            
        }
        
        contigBoundaries = newSeqContigBoundaries;
        
    }
    
    void setSeqGapBoundaries(std::vector<unsigned int> &g) {
        gapBoundaries = g;
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
    
    unsigned int getSeqScaffLen() {
        return inSequence.size();
    }
    
    std::vector<unsigned int> getSeqContigBoundaries() {
        return contigBoundaries;
    }
    
    std::vector<unsigned int> getSeqContigLens() {
        return bedIntervalSizes(contigBoundaries);
    }
    
    unsigned int getContigSum() {
        
        unsigned int contigSum = 0;
        
        for (unsigned int& g : bedIntervalSizes(contigBoundaries))
            contigSum += g;
        
        return contigSum;
    }
    
    unsigned int getContigN() {
        
        return bedIntervalSizes(contigBoundaries).size();
    }
    
    std::vector<unsigned int> getSeqGapBoundaries() {
        return gapBoundaries;
    }
    
    std::vector<unsigned int> getSeqGapLens() {
        
        return bedIntervalSizes(gapBoundaries);
    }
    
    unsigned int getGapSum() {
        
        unsigned int gapSum = 0;
        
        for (auto& g : bedIntervalSizes(gapBoundaries))
            gapSum += g;
        
        return gapSum;
    }
    
    unsigned int getGapN() {
        
        return bedIntervalSizes(gapBoundaries).size();
    }
    
    void setACGT(unsigned int a, unsigned int c, unsigned int g, unsigned int t) {
        
        A = a;
        C = c;
        G = g;
        T = t;
        
    }
    
    void setLowerCount(unsigned int C) {
        
        lowerCount = C;
        
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
    
    double computeGCcontent() {
        
        double GCcontent = (double) (G + C) / (G + C + A + T) * 100;
        
        return GCcontent;
    }
    
};

class InSequences {
    
private:
    std::vector<InSequence> newSeq = std::vector<InSequence>();
    
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
    
    double scaffauN = 0, scaffauNG = 0, contigauN = 0, contigauNG = 0, gapauN = 0;
    
    InSequence inSequence;
    
    unsigned long long int totScaffLen = 0;
    unsigned long long int totContigLen = 0;
    
    unsigned int
    totGapLen = 0,
    totGapN = 0;
    
    unsigned long int totA = 0;
    unsigned long int totC = 0;
    unsigned long int totG = 0;
    unsigned long int totT = 0;
    unsigned long int totLowerCount = 0;
    
public:
    void appendSequence(std::string* seqHeader, std::string* seqComment, std::string* sequence, std::string* sequenceQuality = NULL) {
        
        inSequence.setSeqHeader(seqHeader);
        
        if (seqComment != NULL) {
            
            inSequence.setSeqComment(*seqComment);
            
        }
        
        verbose(verbose_flag, "Header, comment, and sequence read");
        
        verbose(verbose_flag, "Processing scaffold: " + *seqHeader);
        
        inSequence.setInSequence(sequence);
        
        verbose(verbose_flag, "Sequence set");
        
        inSequence.traverseInSequence(sequence);
        
        verbose(verbose_flag, "Sequence traversed");
        
        if (sequenceQuality != NULL) {
        
            inSequence.setInSequenceQuality(sequenceQuality);
        
            verbose(verbose_flag, "Sequence quality set");
            
        }
        
        newSeq.push_back(inSequence);
        
        verbose(verbose_flag, "Sequence added to sequence vector");
        
        increaseTotScaffLen(inSequence.getSeqScaffLen());
        
        verbose(verbose_flag, "Increased total scaffold length");
        
        recordScaffLen(inSequence.getSeqScaffLen());
        
        verbose(verbose_flag, "Recorded length of sequence");
        
        increaseTotContigLen(inSequence.getContigSum());
        
        verbose(verbose_flag, "Increased total contig length");
        
        recordContigLens(inSequence.getSeqContigLens());
        
        verbose(verbose_flag, "Recorded length of contigs in sequence");
        
        recordGapLens(inSequence.getSeqGapLens());
        
        verbose(verbose_flag, "Recorded length of gaps in sequence");
        
        increaseTotGapLen(inSequence.getGapSum());
        
        verbose(verbose_flag, "Increased total gap length");
        
        increaseGapN(inSequence.getGapN());
        
        verbose(verbose_flag, "Increased total number of gaps");
        
        increaseTotACGT(inSequence.getA(), inSequence.getC(), inSequence.getG(), inSequence.getT());
        
        verbose(verbose_flag, "Increased ACGT counts");
        
        increaseTotLowerCount(inSequence.getLowerCount());
        
        verbose(verbose_flag, "Increased count of upper bases");
        
        if(verbose_flag) {std::cout<<"\n";};
        
    }
    
    InSequence getInSequence(unsigned int &idx) {
        
        InSequence inSequence = newSeq[idx];
        return inSequence;
        
    }
    
    void increaseTotScaffLen(unsigned int ScaffLen) {
        
        totScaffLen += ScaffLen;
        
    }
    
    unsigned long long int getTotScaffLen() {
        
        return totScaffLen;
        
    }
    
    void increaseTotContigLen(unsigned int contigLen) {
        
        totContigLen += contigLen;
        
    }
    
    unsigned long long int getTotContigLen() {
        
        return totContigLen;
        
    }
    
    void increaseTotGapLen(unsigned int gapLen) {
        
        totGapLen += gapLen;
        
    }
    
    unsigned int getTotGapLen() {
        
        return totGapLen;
        
    }
    
    void increaseGapN(unsigned int gapN) {
        
        totGapN += gapN;
        
    }
    
    unsigned int getTotGapN() {
        
        return totGapN;
        
    }
    
    void recordScaffLen(unsigned int seqLen) {
        
        scaffLens.push_back(seqLen);
        
    }
    
    void recordContigLens(std::vector <unsigned int> seqLens) {
        
        std::vector <unsigned int> newContigLens;
        
        newContigLens.reserve(contigLens.size() + seqLens.size());
        newContigLens.insert(newContigLens.end(), contigLens.begin(), contigLens.end());
        newContigLens.insert(newContigLens.end(), seqLens.begin(), seqLens.end());
        
        contigLens = newContigLens;
        
    }
    
    void recordGapLens(std::vector <unsigned int> seqLens) {
        
        std::vector <unsigned int> newGapLens;
        
        newGapLens.reserve(gapLens.size() + seqLens.size());
        newGapLens.insert(newGapLens.end(), gapLens.begin(), gapLens.end());
        newGapLens.insert(newGapLens.end(), seqLens.begin(), seqLens.end());
        
        gapLens = newGapLens;
        
    }
    
    void computeScaffNstars(unsigned int gSize) {
        
        sort(scaffLens.begin(), scaffLens.end(), std::greater<unsigned int>());
        
        unsigned long long int scaffSum = 0;
        
        double N = 1, NG = 1;
        
        for(unsigned int i = 0; i < getScaffN(); i++) {
            
            scaffSum += scaffLens[i];
            
            while (scaffSum >= ((double) getTotScaffLen() / 10 * N) && N<= 10) {
                
                scaffNstars[N-1] = scaffLens[i];
                scaffLstars[N-1] = i + 1;
                
                N = N + 1;
                
            }
            
            
            while (gSize > 0 && (scaffSum >= ((double) gSize / 10 * NG)) && NG<= 10) {
                
                scaffNGstars[NG-1] = scaffLens[i];
                scaffLGstars[NG-1] = i + 1;
                
                NG = NG + 1;
                
            }
            
        }
        
    }
    
    void computeContigNstars(unsigned int gSize) {
        
        sort(contigLens.begin(), contigLens.end(), std::greater<unsigned int>());
        
        unsigned long long int contigSum = 0;
        
        short int N = 1, NG = 1;
        
        for(unsigned int i = 0; i < contigLens.size(); i++) {
            
            contigSum += contigLens[i];
            
            while (contigSum >= ((double) getTotContigLen() / 10 * N) && N<= 10) {
                
                contigNstars[N-1] = contigLens[i];
                contigLstars[N-1] = i + 1;
                
                N = N + 1;
                
            }
            
            while (gSize > 0 && (contigSum >= ((double) gSize / 10 * NG)) && NG<= 10) {
                
                contigNGstars[NG-1] = contigLens[i];
                contigLGstars[NG-1] = i + 1;
                
                NG = NG + 1;
                
            }
            
        }
        
    }
    
    void computeGapNstars() {
        
        sort(gapLens.begin(), gapLens.end(), std::greater<unsigned int>());
        
        unsigned long long int gapSum = 0;
        
        short int N = 1;
        
        for(unsigned int i = 0; i < gapLens.size(); i++) {
            
            gapSum += gapLens[i];
            
            while (gapSum >= ((double) totGapLen / 10 * N) && N<= 10) {
                
                gapNstars[N-1] = gapLens[i];
                gapLstars[N-1] = i + 1;
                
                N = N + 1;
                
            }
            
        }
        
    }
    
    void computeScaffauNstar(unsigned int gSize) {
        
        unsigned long long int scaffSum = getTotScaffLen();
        
        for(unsigned int i = 0; i < getScaffN(); i++) {
            
            scaffauN += (double) scaffLens[i] * scaffLens[i] / scaffSum;
            
            if (gSize > 0) {
            
                scaffauNG += (double) scaffLens[i] * scaffLens[i] / gSize;
                
            }
            
        }
        
    }

    void computeContigauNstar(unsigned int gSize) {
        
        unsigned long long int contigSum = getTotContigLen();
        
        for(unsigned int i = 0; i < contigLens.size(); i++) {
            
            contigauN += (double) contigLens[i] * contigLens[i] / contigSum;
            
            if (gSize > 0) {
            
                contigauNG += (double) contigLens[i] * contigLens[i] / gSize;
                
            }
            
        }
        
    }
    
    void computeGapauNstar() {
        
        unsigned long long int gapSum = getTotGapLen();
        
        for(unsigned int i = 0; i < gapLens.size(); i++) {
            
            gapauN += (double) gapLens[i] * gapLens[i] / gapSum;
            
        }
        
    }
    
    unsigned int getScaffN() {
        
        return newSeq.size();
        
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
        
        return scaffauN;
        
    }
    
    double getScaffauNG() {
        
        return scaffauNG;
        
    }
    
    double getContigauN() {
        
        return contigauN;
        
    }
    
    double getContigauNG() {
        
        return contigauNG;
        
    }
    
    double getGapauN() {
        
        return gapauN;
        
    }
    
    unsigned int getContigN() {
        
        return contigLens.size();
        
    }
    
    unsigned int getContigN50() {
        
        return contigNstars[4];
        
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
        
        if (newSeq.size()>0) {
            
            return scaffLens[0];
            
        }else{
            
            return 0;
            
        }
        
    }
    
    double computeAverageScaffLen() {
        
        double AverageScaffLen = (double) totScaffLen/scaffLens.size();
        
        return AverageScaffLen;
        
    }
    
    void increaseTotACGT(unsigned int A, unsigned int C, unsigned int G, unsigned int T) {
        
        totA += A;
        totC += C;
        totG += G;
        totT += T;
        
    }
    
    void increaseTotLowerCount(unsigned int C) {
        
        totLowerCount += C;
        
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
        
        if (newSeq.size()>0) {
            
            GCcontent = (double) (totC + totG) / (totA + totC + totG + totT) * 100;
            
        }else{
            
            GCcontent = 0;
            
        }
        
        return GCcontent;
    }
    
};


#endif /* gfastats-commons_h */
