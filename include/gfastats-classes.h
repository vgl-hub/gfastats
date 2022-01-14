//
//  gfastats-classes.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef gfastatsClasses_h
#define gfastatsClasses_h

//classes
class BedCoordinates {
private:
    std::vector<std::string> fastaHeaders;
    std::vector<unsigned int> cBegin;
    std::vector<unsigned int> cEnd;
    
public:
    
    void pushCoordinates(std::string h, unsigned int b, unsigned int e) {
        
        fastaHeaders.push_back(h);
        cBegin.push_back(b);
        cEnd.push_back(e);
        
    }
    
    bool empty() {
        
        return (fastaHeaders.size()==0) ? true : false;
        
    }
    
    std::vector<std::string> getFastaHeaders() {
        
        return fastaHeaders;
        
    }
    
    std::string getFastaHeader(unsigned int pos) {
        
        return fastaHeaders[pos];
        
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
    std::string fastaHeader;
    std::string fastaComment;
    std::string inSequence;
    std::vector<unsigned int> fastaContigBoundaries;
    std::vector<unsigned int> fastaGapBoundaries;
    unsigned int A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
    
public:
    
    void traverseInSequence(std::string* s) {
        
        unsigned int pos = 0, A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
        bool wasN = false, pushbackGap = false;
        std::vector<unsigned int> caseBoundaries;
        std::vector<unsigned int> fastaGapBoundaries;
        fastaGapBoundaries.reserve(200);
            
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
                            
                            fastaGapBoundaries.push_back(pos);
                            
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
                
                fastaGapBoundaries.push_back(pos);
                pushbackGap = false;
                
            }
            
            pos++;
            
        }
        
        setFastaGapBoundaries(fastaGapBoundaries);
        setFastaContigBoundaries(fastaGapBoundaries);
        setACGT(A, C, G, T);
        setLowerCount(lowerCount);
        
    }
    
    void setFastaHeader(std::string* h) {
        fastaHeader = *h;
    }
    
    void setFastaComment(std::string c) {
        fastaComment = c;
    }
    
    void setInSequence(std::string* s) {
        inSequence = *s;
    }
    
    void setFastaContigBoundaries(std::vector<unsigned int> &fastaGapBoundaries) {
        
        std::vector<unsigned int> newFastaContigBoundaries;
        
        newFastaContigBoundaries.reserve(fastaGapBoundaries.size() + 2);
        
        if (fastaGapBoundaries.size() > 0) {
            
            newFastaContigBoundaries = fastaGapBoundaries;
            
            if (fastaGapBoundaries[0] != 0) {
                
                newFastaContigBoundaries.insert(newFastaContigBoundaries.begin(), 0);
                
            }else{
                
                newFastaContigBoundaries.erase(newFastaContigBoundaries.begin());
                
            }
            
            if (newFastaContigBoundaries[newFastaContigBoundaries.size()-1] != inSequence.size()) {
                
                newFastaContigBoundaries.insert(newFastaContigBoundaries.end(), inSequence.size());
                
            }else{
                
                newFastaContigBoundaries.pop_back();
                
            }
            
        }else{
            
            newFastaContigBoundaries = {0, (unsigned int) inSequence.size()};
            
        }
        
        fastaContigBoundaries = newFastaContigBoundaries;
        
    }
    
    void setFastaGapBoundaries(std::vector<unsigned int> &g) {
        fastaGapBoundaries = g;
    }
    
    std::string getFastaHeader() {
        return fastaHeader;
    }
    
    std::string getFastaComment() {
        return fastaComment;
    }
    
    std::string getInSequence() {
        return inSequence;
    }
    
    unsigned int getFastaScaffLen() {
        return inSequence.size();
    }
    
    std::vector<unsigned int> getFastaContigBoundaries() {
        return fastaContigBoundaries;
    }
    
    std::vector<unsigned int> getFastaContigLens() {
        return bedIntervalSizes(fastaContigBoundaries);
    }
    
    unsigned int getContigSum() {
        
        unsigned int contigSum = 0;
        
        for (unsigned int& g : bedIntervalSizes(fastaContigBoundaries))
            contigSum += g;
        
        return contigSum;
    }
    
    unsigned int getContigN() {
        
        return bedIntervalSizes(fastaContigBoundaries).size();
    }
    
    std::vector<unsigned int> getFastaGapBoundaries() {
        return fastaGapBoundaries;
    }
    
    std::vector<unsigned int> getFastaGapLens() {
        
        return bedIntervalSizes(fastaGapBoundaries);
    }
    
    unsigned int getGapSum() {
        
        unsigned int gapSum = 0;
        
        for (auto& g : bedIntervalSizes(fastaGapBoundaries))
            gapSum += g;
        
        return gapSum;
    }
    
    unsigned int getGapN() {
        
        return bedIntervalSizes(fastaGapBoundaries).size();
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
    std::vector<InSequence> newFasta = std::vector<InSequence>();
    
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
    
    std::string h;
    
public:
    void appendFasta(std::string* h, std::string* c, std::string* s) {
        
        inSequence.setFastaHeader(h);
        
        if (c != NULL) {
            
            inSequence.setFastaComment(*c);
            
        }
        
        verbose(verbose_flag, "Header, comment, and fasta sequence read");
        
        verbose(verbose_flag, "Processing scaffold: " + *h);
        
        inSequence.setInSequence(s);
        
        verbose(verbose_flag, "Fasta sequence set");
        
        inSequence.traverseInSequence(s);
        
        verbose(verbose_flag, "Traversed fasta sequence");
        
        newFasta.push_back(inSequence);
        
        verbose(verbose_flag, "Fasta sequence added to fasta sequence std::vector");
        
        increaseTotScaffLen(inSequence.getFastaScaffLen());
        
        verbose(verbose_flag, "Increased total scaffold length");
        
        recordScaffLen(inSequence.getFastaScaffLen());
        
        verbose(verbose_flag, "Recorded length of fasta sequence");
        
        increaseTotContigLen(inSequence.getContigSum());
        
        verbose(verbose_flag, "Increased total contig length");
        
        recordContigLens(inSequence.getFastaContigLens());
        
        verbose(verbose_flag, "Recorded length of contigs in fasta sequence");
        
        recordGapLens(inSequence.getFastaGapLens());
        
        verbose(verbose_flag, "Recorded length of gaps in fasta sequence");
        
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
        
        InSequence fastaSequence = newFasta[idx];
        return fastaSequence;
        
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
        
        return newFasta.size();
        
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
        
        if (newFasta.size()>0) {
            
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
        
        if (newFasta.size()>0) {
            
            GCcontent = (double) (totC + totG) / (totA + totC + totG + totT) * 100;
            
        }else{
            
            GCcontent = 0;
            
        }
        
        return GCcontent;
    }
    
};

class InFile {
    
    std::string h;
    char* c;
    std::vector<std::string> bedIncludeListHeaders;
    std::vector<std::string> bedExcludeListHeaders;
    unsigned int pos = 0, cBegin = 0, cEnd = 0, offset = 0, prevCEnd = 0;
    
public:
    
    InSequences readFiles(std::string &iFastaFileArg, std::string &iBedIncludeFileArg, std::string &iBedExcludeFileArg, BedCoordinates &bedIncludeList, bool isPipe, char &pipeType) {
        
        std::string newLine, fastaHeader, fastaComment, fastaSequence, line, h;
        std::unique_ptr<std::istream> stream;
        
        unsigned int idx = 0, b = 0, e = 0;
        
        if (!iBedIncludeFileArg.empty() || (isPipe && (pipeType == 'i'))) {
            
            if (isPipe && (pipeType == 'i')) {
                
                std::istream &in = std::cin;
                stream = make_unique<std::istream>(in.rdbuf());
                
            }else{
                
                stream = make_unique<std::ifstream>(std::ifstream(iBedIncludeFileArg));
                
            }
            
            while (getline(*stream, line)) {
                
                std::istringstream iss(line);
                iss >> h >> b >> e;
                
                bedIncludeList.pushCoordinates(h, b, e);
                b = 0, e = 0;
                
            }
            
        }
        
        BedCoordinates bedExcludeList;
        
        if (!iBedExcludeFileArg.empty() || (isPipe && (pipeType == 'e'))) {
            
            if (isPipe && (pipeType == 'e')) {
                
                std::istream &in = std::cin;
                stream = make_unique<std::istream>(in.rdbuf());
                
            }else{
                
                stream = make_unique<std::ifstream>(std::ifstream(iBedExcludeFileArg));
                
            }
            
            while (getline(*stream, line)) {
                
                std::istringstream iss(line);
                iss >> h >> b >> e;
                
                bedExcludeList.pushCoordinates(h, b, e);
                b = 0, e = 0;
                
            }
            
        }
        
        InSequences Fasta;
        
        std::string firstLine;
        char firstChar;
        
        if (determineGzip(iFastaFileArg)) {
            
            std::string data;
            
            data = loadGzip(iFastaFileArg);
            
            stream = make_unique<std::istringstream>(std::istringstream(data));
            
        } else if (isPipe && (pipeType == 'f')) {
            
            std::istream &in = std::cin;
            stream = make_unique<std::istream>(in.rdbuf());
            
        } else {
            
            stream = make_unique<std::ifstream>(std::ifstream(iFastaFileArg));
            
        }
        
        if (stream) {
            
            getline(*stream, newLine);
            firstLine = newLine;
            firstChar = newLine[0];
            
            if (!isPipe || pipeType != 'f') {
                
                stream->clear();
                stream->seekg(0, stream->beg);
                
            }
            
            switch (firstChar) {
                    
                case '>': {
                    
                    if (isPipe && pipeType == 'f') {
                        
                        parseFasta(firstLine, Fasta, fastaHeader, fastaComment, fastaSequence, idx, bedIncludeList, bedExcludeList);
                        
                    }
                    
                    while (getline(*stream, newLine)) {
                        
                        parseFasta(newLine, Fasta, fastaHeader, fastaComment, fastaSequence, idx, bedIncludeList, bedExcludeList);
                        
                    }
                    
                    includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, bedIncludeList, bedExcludeList);
                    
                    break;
                }
                case '@': {
                    
                    if (isPipe && pipeType == 'f') {
                        
                        firstLine.erase(0, 1);
                        
                        h = std::string(strtok(strdup(firstLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        fastaHeader = h;
                        
                        if (c != NULL) {
                            
                            fastaComment = std::string(c);
                            
                        }
                        
                        getline(*stream, newLine);
                        
                        fastaSequence = newLine;
                        
                        getline(*stream, newLine);
                        getline(*stream, newLine);
                        
                        includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, bedIncludeList, bedExcludeList);
                        
                    }
                    
                    while (getline(*stream, newLine)) {
                        
                        newLine.erase(0, 1);
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        fastaHeader = h;
                        
                        if (c != NULL) {
                            
                            fastaComment = std::string(c);
                            
                        }
                        
                        getline(*stream, newLine);
                        
                        fastaSequence = newLine;
                        
                        getline(*stream, newLine);
                        getline(*stream, newLine);
                        
                        includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, bedIncludeList, bedExcludeList);
                        
                    }
                    
                    break;
                    
                }
                default: {
                    
                    std::string h_col1, h_col2, s;
                    char* version;
                    
                    h_col1 = std::string(strtok(strdup(firstLine.c_str()),"\t")); //process first line
                    
                    if(h_col1 == "H"){
                        
                        h_col2 = strtok(NULL,""); //read header col2
                        std::string(strtok(strdup(h_col2.c_str()),":")); //process version tag
                        strtok(NULL,":");
                        version = strtok(NULL,"");
                        
                        if (version != NULL) {
                            
                            verbose(verbose_flag, "GFA version: " + std::string(version));
                            
                        }else{
                            
                            printf("Cannot recognize GFA version");
                            printf("Offending line: %s", firstLine.c_str());
                            exit(1);
                            
                        }
                        
                    }
                    
                    while (getline(*stream, newLine)) {
                        
                        switch (newLine[0]) {
                                
                            case 'S': {
                                
                                std::string(strtok(strdup(newLine.c_str()),"\t")); //process first line
                                h = strtok(NULL,"\t");
                                
                                fastaHeader = h;
                                
                                s = strtok(NULL,"\t");
                                fastaSequence = s;
                                
                                includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, bedIncludeList, bedExcludeList);
                                
                                break;
                                
                            }
                            default: {
                                
                                break;
                                
                            }
                                
                        }
                        
                    }
                    
                    break;
                    
                }
            }
            
        }else{
            
            printf("Stream not successful: %s", iFastaFileArg.c_str());
            
        }
        
        return Fasta;
        
    }
    
    void parseFasta(std::string &newLine, InSequences &Fasta, std::string &fastaHeader, std::string &fastaComment, std::string &inSequence, unsigned int &idx, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList) {
        
        switch (newLine[0]) {
                
            case '>': {
                
                if (idx> 0) {
                    
                    includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &inSequence, bedIncludeList, bedExcludeList);
                    
                    inSequence = "";
                    
                }
                
                newLine.erase(0, 1);
                
                h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                c = strtok(NULL,""); //read comment
                
                fastaHeader = h;
                
                if (c != NULL) {
                    
                    fastaComment = std::string(c);
                    
                }
                
                idx++;
                
                break;
            }
            case '\n':{
                
                break;
            }
            case ' ':{
                
                break;
            }
            default: {
                
                inSequence.append(newLine);
                
            }
                
        }
        
    }
    
    void includeExcludeAppend(InSequences* Fasta, std::string* fastaHeader, std::string* fastaComment, std::string* fastaSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList) {
 
        bedIncludeListHeaders = bedIncludeList.getFastaHeaders();
        bedExcludeListHeaders = bedExcludeList.getFastaHeaders();
        bool outFasta;
        
        if   (bedIncludeList.empty() &&
              bedExcludeList.empty()) {
            
            Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
            
        }else if
            (!bedIncludeList.empty() &&
             bedExcludeList.empty()) {
            
            offset = 0, prevCEnd = 0;
            outFasta = false;
            
            auto it = begin(bedIncludeListHeaders);
            
            while (it != end(bedIncludeListHeaders)) {
                
                it = std::find(it, bedIncludeListHeaders.end(), *fastaHeader);
                
                if (it == end(bedIncludeListHeaders) || bedIncludeList.getFastaHeader(pos) != *fastaHeader) {
                    
                    break;
                    
                }
                
                outFasta = true;

                cBegin = bedIncludeList.getcBegin(pos);
                cEnd = bedIncludeList.getcEnd(pos);
                
                if (!(cBegin == 0 && cEnd == 0)) {
                    
                    fastaSequence->erase(offset, cBegin-prevCEnd);
                    offset += cEnd-cBegin;
                    prevCEnd = cEnd;
                    
                }
              
                ++it;
                pos++;
                
            }
                
            if (outFasta && fastaSequence->size()>0) {
                
                if (offset>0) {
                
                    fastaSequence->erase(offset, fastaSequence->size()-offset);
                    
                }
                
                Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
            
            }else {
                
                verbose(verbose_flag, "Scaffold entirely removed as a result of include: " + *fastaHeader);
                
            }
                
        }else if
            (bedIncludeList.empty() &&
             !bedExcludeList.empty()) {
                
            offset = 0;
            outFasta = true;
            
            auto it = begin(bedExcludeListHeaders);
            
            while (it != end(bedExcludeListHeaders)) {
                
                it = std::find(it, bedExcludeListHeaders.end(), *fastaHeader);
                
                if (it == end(bedExcludeListHeaders) || bedExcludeList.getFastaHeader(pos) != *fastaHeader) {
                    
                    break;
                    
                }

                cBegin = bedExcludeList.getcBegin(pos);
                cEnd = bedExcludeList.getcEnd(pos);
                
                if (!(cBegin == 0 && cEnd == 0)) {
                    
                    fastaSequence->erase(cBegin-offset, cEnd-cBegin);
                    offset += cEnd-cBegin;
                    
                }else{
                    
                    outFasta = false;
                    
                }
              
                ++it;
                pos++;
                
            }
                
            if (outFasta && fastaSequence->size()>0) {
            
                Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
            
            }else {
                
                verbose(verbose_flag, "Scaffold entirely removed as a result of exclude: " + *fastaHeader);
                
            }
                    
        }else if
                (!bedIncludeList.empty() &&
                 !bedExcludeList.empty() &&
                 std::find(bedIncludeListHeaders.begin(), bedIncludeListHeaders.end(), *fastaHeader) != bedIncludeListHeaders.end() &&
                 std::find(bedExcludeListHeaders.begin(), bedExcludeListHeaders.end(), *fastaHeader) == bedExcludeListHeaders.end()) {
                    
                    Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
                    
        }
        
    }
    
};


#endif /* gfastats-classes_h */
