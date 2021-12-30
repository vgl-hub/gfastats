//
//  gfastats-classes.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef gfastatsClasses_h
#define gfastatsClasses_h

//classes
class FastaSequence {
private:
    std::string fastaHeader;
    std::string fastaSequence;
    std::string fastaComment;
    std::vector<unsigned int> fastaContigBoundaries;
    std::vector<unsigned int> fastaGapBoundaries;
    unsigned int A = 0, C = 0, G = 0, T = 0;
    
public:
    
    void TraverseFastaSequence(std::string &s) {
        
        unsigned int pos = 0, A = 0, C = 0, G = 0, T = 0;
        bool wasN = false, pushbackGap = false;
        std::vector<unsigned int> fastaGapBoundaries;
        fastaGapBoundaries.reserve(200);
        
        for (char base : s) {
            
            switch (base) {
                    
                case 'N':
                case 'n':
                case 'X': {
                    
                    if (!wasN) { // gap start
                        
                        pushbackGap = true;
                        
                    }
                    
                    if(pos == (s.length()-1)) { // end of scaffold
                        
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
    
    void setFastaHeader(std::string &h) {
        fastaHeader = h;
    }
    
    void setFastaComment(std::string c) {
        fastaComment = c;
    }
    
    void setFastaSequence(std::string &s) {
        fastaSequence = s;
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
    
    unsigned int getFastaSeqLen() {
        return fastaSequence.size();
    }
    
    std::vector<unsigned int> getFastaContigs() {
        return fastaContigBoundaries;
    }
    
    std::vector<unsigned int> getFastaContigLens() {
        
        return bedIntervalSizes(fastaContigBoundaries);
    }
    
    std::vector<unsigned int> getFastaGaps() {
        return fastaGapBoundaries;
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

class FastaSequences {
    
private:
    std::vector<FastaSequence> newFasta = std::vector<FastaSequence>();
    
    std::vector<unsigned int> scaffLens;
    
    std::vector<unsigned int> contigLens;
    
    FastaSequence fastaSequence;
    
    unsigned long long int totScaffLen = 0;
    
    unsigned int
    totGapLen = 0,
    gapN = 0,
    scaffN50 = 0, scaffNG50 = 0,
    scaffL50 = 0, scaffLG50 = 0,
    contigN50 = 0, contigNG50 = 0,
    contigL50 = 0, contigLG50 = 0;
    
    unsigned long int totA = 0;
    unsigned long int totC = 0;
    unsigned long int totG = 0;
    unsigned long int totT = 0;
    
    std::string h;
    char *c;
    
public:
    void appendFasta(std::string hg, std::string s) {
        
        h = std::string(strtok(strdup(hg.c_str())," ")); //process header line
        h.erase(0, 1);
        fastaSequence.setFastaHeader(h);
        
        c = strtok(NULL,""); //process comment line
        
        if (c != NULL) {
            
            fastaSequence.setFastaComment(std::string(c));
            
        }
        
        verbose(verbose_flag, "Header, comment, and fasta sequence read");
        
        verbose(verbose_flag, "Processing scaffold: "+h);
        
        fastaSequence.setFastaSequence(s);
        
        verbose(verbose_flag, "Fasta sequence set");
        
        fastaSequence.TraverseFastaSequence(s);
        
        verbose(verbose_flag, "Traversed fasta sequence");
        
        newFasta.push_back(fastaSequence);
        
        verbose(verbose_flag, "Fasta sequence added to fasta sequence std::vector");
        
        increaseTotScaffLen(fastaSequence.getFastaSeqLen());
        
        verbose(verbose_flag, "Increased total scaffold length");
        
        recordScaffLen(fastaSequence.getFastaSeqLen());
        
        verbose(verbose_flag, "Recorded length of fasta sequence");
        
        recordContigLens(fastaSequence.getFastaContigLens());
        
        verbose(verbose_flag, "Recorded length of contigs in fasta sequence");
        
        increaseTotGapLen(fastaSequence.getGapSum());
        
        verbose(verbose_flag, "Increased total gap length");
        
        increaseGapN(fastaSequence.getGapN());
        
        verbose(verbose_flag, "Increased total number of gaps");

        increaseTotACGT(fastaSequence.getA(), fastaSequence.getC(), fastaSequence.getG(), fastaSequence.getT());
        
        verbose(verbose_flag, "Increased ACGT counts");
        
        if(verbose_flag) {std::cout<<"\n";};
        
    }
    
    FastaSequence getFastaSequences(unsigned int idx) {
        
        FastaSequence fastaSequence = newFasta[idx];
        return fastaSequence;
        
    }
    
    void increaseTotScaffLen(unsigned int ScaffLen) {
        
        totScaffLen += ScaffLen;
        
    }
    
    unsigned long long int getTotScaffLen() {
        
        return totScaffLen;
        
    }
    
    void increaseTotGapLen(unsigned int GapLen) {
        
        totGapLen += GapLen;
        
    }
    
    unsigned int getTotGapLen() {
        
        return totGapLen;
        
    }
    
    void increaseGapN(unsigned int GapN) {
        
        gapN += GapN;
        
    }
    
    unsigned int getTotGapN() {
        
        return gapN;
        
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
    
    void computeScaffN50(unsigned int gSize) {
        
        sort(scaffLens.begin(), scaffLens.end(), std::greater<unsigned int>());
        
        unsigned long long int scaffSum = 0;
        
        for(unsigned int i = 0; i < getScaffN(); i++) {
            
            scaffSum += scaffLens[i];
            
            if (scaffSum >= getTotScaffLen() / 2 && scaffN50 < scaffLens[i]) {
                
                scaffN50 = scaffLens[i];
                scaffL50 = i;
                
            }
            
            if (gSize > 0 && scaffSum >= gSize / 2 && scaffNG50 < scaffLens[i]) {
                
                scaffNG50 = scaffLens[i];
                scaffLG50 = i;
            }
            
            if (scaffN50 >= scaffLens[i] && scaffNG50 >= scaffLens[i]) {
                
                break;
                
            }
            
        }
        
    }
    
    void computeContigN50(unsigned int gSize) {
        
        sort(contigLens.begin(), contigLens.end(), std::greater<unsigned int>());
        
        unsigned long long int contigSum = 0;
        
        for(unsigned int i = 0; i < contigLens.size(); i++) {
            
            contigSum += contigLens[i];
            
            if (contigSum >= getTotScaffLen() / 2 && contigN50 < contigLens[i]) {
                
                contigN50 = contigLens[i];
                contigL50 = i;
                
            }
            
            if (gSize > 0 && contigSum >= gSize / 2 && contigNG50 < contigLens[i]) {
                
                contigNG50 = contigLens[i];
                contigLG50 = i;
            }
            
            if (contigN50 >= contigLens[i] && contigNG50 >= contigLens[i]) {
                
                break;
                
            }
            
        }
        
    }
    
    unsigned int getScaffN() {
        
        return newFasta.size();
        
    }
    
    unsigned int getScaffN50(unsigned long long int gSize) {
        
        computeScaffN50(gSize);
        
        return scaffN50;
        
    }
    
    unsigned int getScaffNG50() {
        
        return scaffNG50;
        
    }
    
    unsigned int getScaffL50() {
        
        return scaffL50;
        
    }
    
    unsigned int getScaffLG50() {
        
        return scaffLG50;
        
    }
    
    unsigned int getContigN() {
        
        return contigLens.size();
        
    }
    
    unsigned int getContigN50(unsigned long long int gSize) {
        
        computeContigN50(gSize);
        
        return contigN50;
        
    }

    unsigned int getContigNG50() {
        
        return contigNG50;
        
    }

    unsigned int getContigL50() {
        
        return contigL50;
        
    }
    
    unsigned int getContigLG50() {
        
        return contigLG50;
        
    }
    
    unsigned int getLargestScaffold() {
        
        return scaffLens[0];
        
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
        
        double GCcontent = (double) (totC + totG) / (totA + totC + totG + totT) * 100;
    
        return GCcontent;
    }
    
};

class FastaFile {
    
public:
    
    void ParseFasta(std::string newLine, FastaSequences &Fasta, std::string &fastaHeader, std::string &fastaSequence, unsigned int &idx) {
        
        switch (newLine[0]) {
                
            case '>': {
                
                if (idx> 0) {
                    
                    Fasta.appendFasta(fastaHeader,fastaSequence);
                    fastaSequence = "";
                    
                }
                
                fastaHeader = newLine;
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
    
    
    FastaSequences Read(std::string iFileArg) {
        
        std::string newLine, fastaHeader, fastaSequence;
        unsigned int idx = 0;
        
        FastaSequences Fasta;
        
        std::ifstream stream(iFileArg);
        
        unsigned char buffer[2];
        stream.read((char*)(&buffer[0]), 2) ;
        
        
        stream.clear();
        stream.seekg(0, stream.beg);
        
        if (buffer[0] == 0x1f && (buffer[1] == 0x8b)) {
            
            stream.close();
            
            std::string fileData;
            if (!loadBinaryFile(iFileArg, fileData)) {
                printf("Error loading input file.");
            }
            
            std::string data;
            if (!gzipInflate(fileData, data)) {
                printf("Error decompressing input file.");
                
            }
            
            std::stringstream gzstream(data);
            
            while (getline (gzstream, newLine)) {
                
                if (gzstream) {
                    
                    ParseFasta(newLine, Fasta, fastaHeader, fastaSequence, idx);
                    
                }
                
                else {
                    
                    std::cout << "Gzip stream not successful.";
                    
                }
            }
            
        } else {
            
            while (getline (stream, newLine)) {
                
                ParseFasta(newLine, Fasta, fastaHeader, fastaSequence, idx);
                
            }
            
            stream.close();
            
        }
        
        Fasta.appendFasta(fastaHeader,fastaSequence);
        
        return Fasta;
        
    }
};


#endif /* gfastats-classes_h */
