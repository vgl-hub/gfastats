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
        
        for (char &base : s) {
            
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
    
    unsigned int getFastaScaffLens() {
        return fastaSequence.size();
    }
    
    std::vector<unsigned int> getFastaContigs() {
        return fastaContigBoundaries;
    }
    
    std::vector<unsigned int> getFastaContigLens() {
        
        return bedIntervalSizes(fastaContigBoundaries);
    }
    
    unsigned int getContigSum() {
        
        unsigned int contigSum = 0;
        
        for (auto& g : bedIntervalSizes(fastaContigBoundaries))
            contigSum += g;
        
        return contigSum;
    }

    std::vector<unsigned int> getFastaGaps() {
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

class FastaSequences {
    
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
    gapN = 0;
    
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
        
        verbose(verbose_flag, "Processing scaffold: " + h);
        
        fastaSequence.setFastaSequence(s);
        
        verbose(verbose_flag, "Fasta sequence set");
        
        fastaSequence.TraverseFastaSequence(s);
        
        verbose(verbose_flag, "Traversed fasta sequence");
        
        newFasta.push_back(fastaSequence);
        
        verbose(verbose_flag, "Fasta sequence added to fasta sequence std::vector");
        
        increaseTotScaffLen(fastaSequence.getFastaScaffLens());
        
        verbose(verbose_flag, "Increased total scaffold length");
        
        recordScaffLen(fastaSequence.getFastaScaffLens());
        
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
    
    void increaseTotContigLen(unsigned int ContigLen) {
        
        totContigLen += ContigLen;
        
    }
    
    unsigned long long int getTotContigLen() {
        
        return totContigLen;
        
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
        
        std::cout<<"scaffLens: "<<getScaffN()<<std::endl;
        
        for(unsigned int i = 0; i < getScaffN(); i++) {
            
            scaffSum += scaffLens[i];

            if (scaffSum >= ((double) getTotScaffLen() / 10 * N) && N<= 10) {
                
                scaffNstars[N-1] = scaffLens[i];
                scaffLstars[N-1] = i + 1;
                
                N = N + 1;
                
            }
            
            
            if (gSize > 0 && (scaffSum >= ((double) gSize / 10 * NG)) && NG<= 10) {
                
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
            
            if (contigSum >= ((double) getTotContigLen() / 10 * N) && N<= 10) {
                
                contigNstars[N-1] = contigLens[i];
                contigLstars[N-1] = i + 1;
                
                N = N + 1;
                
            }
            
            if (gSize > 0 && (contigSum >= ((double) gSize / 10 * NG)) && NG<= 10) {
                
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
