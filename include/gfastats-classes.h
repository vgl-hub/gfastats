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

class FastaSequence {
private:
    std::string fastaHeader;
    std::string fastaSequence;
    std::string fastaComment;
    std::vector<unsigned int> fastaContigBoundaries;
    std::vector<unsigned int> fastaGapBoundaries;
    unsigned int A = 0, C = 0, G = 0, T = 0;
    
public:
    
    void TraverseFastaSequence(std::string* s) {
        
        unsigned int pos = 0, A = 0, C = 0, G = 0, T = 0;
        bool wasN = false, pushbackGap = false;
        std::vector<unsigned int> fastaGapBoundaries;
        fastaGapBoundaries.reserve(200);
        
        for (char &base : *s) {
            
            switch (base) {
                    
                case 'N':
                case 'n':
                case 'X': {
                    
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
        
    }
    
    void setFastaHeader(std::string* h) {
        fastaHeader = *h;
    }
    
    void setFastaComment(std::string c) {
        fastaComment = c;
    }
    
    void setFastaSequence(std::string* s) {
        fastaSequence = *s;
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
            
            if (newFastaContigBoundaries[newFastaContigBoundaries.size()-1] != fastaSequence.size()) {
                
                newFastaContigBoundaries.insert(newFastaContigBoundaries.end(), fastaSequence.size());
                
            }else{
                
                newFastaContigBoundaries.pop_back();
                
            }
            
        }else{
            
            newFastaContigBoundaries = {0, (unsigned int) fastaSequence.size()};
            
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
    
    std::string getFastaSequence() {
        return fastaSequence;
    }
    
    unsigned int getFastaScaffLen() {
        return fastaSequence.size();
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
    
    double computeGCcontent() {
        
        double GCcontent = (double) (G + C) / (G + C + A + T) * 100;
        
        return GCcontent;
    }
    
};

class iSequences {
    
private:
    std::vector<FastaSequence> newFasta = std::vector<FastaSequence>();
    
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
    
    FastaSequence fastaSequence;
    
    unsigned long long int totScaffLen = 0;
    unsigned long long int totContigLen = 0;
    
    unsigned int
    totGapLen = 0,
    totGapN = 0;
    
    unsigned long int totA = 0;
    unsigned long int totC = 0;
    unsigned long int totG = 0;
    unsigned long int totT = 0;
    
    std::string h;
    
public:
    void appendFasta(std::string* h, std::string* c, std::string* s) {
        
        fastaSequence.setFastaHeader(h);
        
        if (c != NULL) {
            
            fastaSequence.setFastaComment(*c);
            
        }
        
        verbose(verbose_flag, "Header, comment, and fasta sequence read");
        
        verbose(verbose_flag, "Processing scaffold: " + *h);
        
        fastaSequence.setFastaSequence(s);
        
        verbose(verbose_flag, "Fasta sequence set");
        
        fastaSequence.TraverseFastaSequence(s);
        
        verbose(verbose_flag, "Traversed fasta sequence");
        
        newFasta.push_back(fastaSequence);
        
        verbose(verbose_flag, "Fasta sequence added to fasta sequence std::vector");
        
        increaseTotScaffLen(fastaSequence.getFastaScaffLen());
        
        verbose(verbose_flag, "Increased total scaffold length");
        
        recordScaffLen(fastaSequence.getFastaScaffLen());
        
        verbose(verbose_flag, "Recorded length of fasta sequence");
        
        increaseTotContigLen(fastaSequence.getContigSum());
        
        verbose(verbose_flag, "Increased total contig length");
        
        recordContigLens(fastaSequence.getFastaContigLens());
        
        verbose(verbose_flag, "Recorded length of contigs in fasta sequence");
        
        recordGapLens(fastaSequence.getFastaGapLens());
        
        verbose(verbose_flag, "Recorded length of gaps in fasta sequence");
        
        increaseTotGapLen(fastaSequence.getGapSum());
        
        verbose(verbose_flag, "Increased total gap length");
        
        increaseGapN(fastaSequence.getGapN());
        
        verbose(verbose_flag, "Increased total number of gaps");
        
        increaseTotACGT(fastaSequence.getA(), fastaSequence.getC(), fastaSequence.getG(), fastaSequence.getT());
        
        verbose(verbose_flag, "Increased ACGT counts");
        
        if(verbose_flag) {std::cout<<"\n";};
        
    }
    
    FastaSequence getISequence(unsigned int &idx) {
        
        FastaSequence fastaSequence = newFasta[idx];
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
    
    void computeGapNstars(unsigned int gSize) {
        
        sort(gapLens.begin(), gapLens.end(), std::greater<unsigned int>());
        
        unsigned long long int gapSum = 0;
        
        short int N = 1;
        
        for(unsigned int i = 0; i < gapLens.size(); i++) {
            
            gapSum += gapLens[i];
            
            if (gapSum >= ((double) totGapLen / 10 * N) && N<= 10) {
                
                gapNstars[N-1] = gapLens[i];
                gapLstars[N-1] = i + 1;
                
                N = N + 1;
                
            }
            
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

class iFile {
    
    std::string h;
    char* c;
    std::vector<std::string> headerIncludeListCpy;
    std::vector<std::string> bedExcludeListHeaders;
    unsigned int pos = 0, cBegin = 0, cEnd = 0, offset = 0;
    
public:
    
    iSequences readFiles(std::string &iFastaFileArg, std::string &iBedIncludeFileArg, std::string &iBedExcludeFileArg, BedCoordinates &headerIncludeList, bool isPipe, char &pipeType) {
        
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
                
                headerIncludeList.pushCoordinates(h, b, e);
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
        
        iSequences Fasta;
        
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
                        
                        parseFasta(firstLine, Fasta, fastaHeader, fastaComment, fastaSequence, idx, headerIncludeList, bedExcludeList);
                        
                    }
                    
                    while (getline(*stream, newLine)) {
                        
                        parseFasta(newLine, Fasta, fastaHeader, fastaComment, fastaSequence, idx, headerIncludeList, bedExcludeList);
                        
                    }
                    
                    includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, headerIncludeList, bedExcludeList);
                    
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
                        
                        includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, headerIncludeList, bedExcludeList);
                        
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
                        
                        includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, headerIncludeList, bedExcludeList);
                        
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
                                
                                includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, headerIncludeList, bedExcludeList);
                                
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
    
    void parseFasta(std::string &newLine, iSequences &Fasta, std::string &fastaHeader, std::string &fastaComment, std::string &fastaSequence, unsigned int &idx, BedCoordinates headerIncludeList, BedCoordinates bedExcludeList) {
        
        switch (newLine[0]) {
                
            case '>': {
                
                if (idx> 0) {
                    
                    includeExcludeAppend(&Fasta, &fastaHeader, &fastaComment, &fastaSequence, headerIncludeList, bedExcludeList);
                    
                    fastaSequence = "";
                    
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
                
                fastaSequence.append(newLine);
                
            }
                
        }
        
    }
    
    void includeExcludeAppend(iSequences* Fasta, std::string* fastaHeader, std::string* fastaComment, std::string* fastaSequence, BedCoordinates headerIncludeList, BedCoordinates bedExcludeList) {
        
        headerIncludeListCpy = headerIncludeList.getFastaHeaders();
        bedExcludeListHeaders = bedExcludeList.getFastaHeaders();
        
        if   (headerIncludeList.empty() &&
              bedExcludeList.empty()) {
            
            Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
            
        }else if
            (!headerIncludeList.empty() &&
             bedExcludeList.empty()) {
                
                auto it = std::find(headerIncludeListCpy.begin(), headerIncludeListCpy.end(), *fastaHeader);
                
                if (it != headerIncludeListCpy.end()) {
                    
                    pos = it - headerIncludeListCpy.begin();
                    
                    cBegin = headerIncludeList.getcBegin(pos);
                    cEnd = headerIncludeList.getcEnd(pos);
                    
                    std::cout<<*fastaSequence<<std::endl;
                    
                    if (!(cBegin == 0 && cEnd == 0)) {
                        
                        std::cout<<cBegin<<std::endl;
                        std::cout<<cEnd<<std::endl;
                        fastaSequence->erase(cEnd, fastaSequence->size()-1);
                        fastaSequence->erase(0, cBegin);
                        
                    }
                    
                    std::cout<<*fastaSequence<<std::endl;
                    
                    pos = 0;
                    Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
                    
                }
                
        }else if
            (headerIncludeList.empty() &&
             !bedExcludeList.empty()) {
                
            offset = 0;
            
            auto it = begin(bedExcludeListHeaders);
            
            while (it != end(bedExcludeListHeaders)) {
                
                it = std::find(it, bedExcludeListHeaders.end(), *fastaHeader);
                
                if (it == end(bedExcludeListHeaders) || bedExcludeList.getFastaHeader(pos) != *fastaHeader) {
                    
                    break;
                    
                }

                cBegin = bedExcludeList.getcBegin(pos);
                cEnd = bedExcludeList.getcEnd(pos);
                
                if (!(cBegin == 0 && cEnd == 0)) {
                    
                    fastaSequence->erase(cBegin-offset, cEnd-cBegin-offset);
                    offset += cEnd-cBegin;
                    
                }
              
                ++it;
                pos++;
                
            }
            
            Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
                    
        }else if
            (headerIncludeList.empty() &&
             !bedExcludeList.empty() &&
             std::find(bedExcludeListHeaders.begin(), bedExcludeListHeaders.end(), *fastaHeader) == bedExcludeListHeaders.end()) {
                
                Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
                
        }else if
                (!headerIncludeList.empty() &&
                 !bedExcludeList.empty() &&
                 std::find(headerIncludeListCpy.begin(), headerIncludeListCpy.end(), *fastaHeader) != headerIncludeListCpy.end() &&
                 std::find(bedExcludeListHeaders.begin(), bedExcludeListHeaders.end(), *fastaHeader) == bedExcludeListHeaders.end()) {
                    
                    Fasta->appendFasta(fastaHeader, fastaComment, fastaSequence);
                    
        }
        
    }
    
};


#endif /* gfastats-classes_h */
