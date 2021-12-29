#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <unistd.h>
#include <getopt.h>

#include <vector>
#include <algorithm>

#include <chrono>

#include "zlib.h"

using namespace std;

//global
static auto start = chrono::high_resolution_clock::now();

//functions
bool isN(char base){
    return (base == 'N' || base == 'n' || base == 'X');
}

double elapsedTime(){
    
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    start = chrono::high_resolution_clock::now();
    
    return elapsed.count();
    
}

void verbose(short int verbose_flag, string msg) {
    
    if(verbose_flag) {
        
        cout << msg<< " (done in " << elapsedTime() << " s).\n";
        
        elapsedTime();
        
    }
}

vector<unsigned int> bedIntervalSizes(vector<unsigned int> intervalVec){
    
    vector<unsigned int> intervalVecLens;
    intervalVecLens.reserve(200);
    
    if (!intervalVec.empty()){
        vector<unsigned int>::const_iterator end = intervalVec.cend();
        
        for (vector<unsigned int>::const_iterator it = intervalVec.cbegin(); it != end;) {
            
            intervalVecLens.push_back(*(it+1) - *it);
            
            it = it + 2;
            
        }
        
    }
    
    return intervalVecLens;
    
}

static short int tabular_flag;
string output(string output){
    
    if (tabular_flag) {
        
        output = output + "\t";
        
    }else{
        
        output = output + ": ";
        
    }
    
    return output;
}

//zlib
bool gzipInflate(const std::string& compressedBytes, std::string& uncompressedBytes) {
    if (compressedBytes.size() == 0) {
        uncompressedBytes = compressedBytes ;
        return true ;
    }
    
    uncompressedBytes.clear() ;
    
    unsigned full_length = compressedBytes.size() ;
    unsigned half_length = compressedBytes.size() / 2;
    
    unsigned uncompLength = full_length ;
    char* uncomp = (char*) calloc(sizeof(char), uncompLength);
    
    z_stream strm;
    strm.next_in = (Bytef *) compressedBytes.c_str();
    strm.avail_in = compressedBytes.size() ;
    strm.total_out = 0;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    
    bool done = false ;
    
    if (inflateInit2(&strm, (16+MAX_WBITS)) != Z_OK) {
        free(uncomp);
        return false;
    }
    
    while (!done) {
        // If our output buffer is too small
        if (strm.total_out >= uncompLength) {
            // Increase size of output buffer
            char* uncomp2 = (char*) calloc(sizeof(char), uncompLength + half_length);
            memcpy(uncomp2, uncomp, uncompLength);
            uncompLength += half_length ;
            free(uncomp);
            uncomp = uncomp2 ;
        }
        
        strm.next_out = (Bytef *) (uncomp + strm.total_out);
        strm.avail_out = uncompLength - strm.total_out;
        
        // Inflate another chunk.
        int err = inflate (&strm, Z_SYNC_FLUSH);
        if (err == Z_STREAM_END) done = true;
        else if (err != Z_OK)  {
            break;
        }
    }
    
    if (inflateEnd (&strm) != Z_OK) {
        free(uncomp);
        return false;
    }
    
    for (size_t i=0; i<strm.total_out; ++i) {
        uncompressedBytes += uncomp[i];
    }
    free(uncomp);
    return true ;
}

/* Reads a file into memory. */
bool loadBinaryFile(const std::string& filename, std::string& contents) {
    // Open the gzip file in binary mode
    FILE* f = fopen(filename.c_str(), "rb");
    if (f == NULL)
        return false ;
    
    // Clear existing bytes in output vector
    contents.clear();
    
    // Read all the bytes in the file
    int c = fgetc(f);
    while (c != EOF) {
        contents +=  (char) c ;
        c = fgetc(f);
    }
    fclose (f);
    
    return true ;
}

//classes
class FastaSequence {
private:
    string fastaHeader;
    string fastaSequence;
    string fastaComment;
    vector<unsigned int> fastaContigBoundaries;
    vector<unsigned int> fastaGapBoundaries;
    
public:
    
    void TraverseFastaSequence(string s) {
        
        unsigned int pos = 0;
        bool wasN = false, pushbackGap = false;
        vector<unsigned int> fastaGapBoundaries;
        fastaGapBoundaries.reserve(200);
        
        for (char base : s) {
            
            if (isN(base) && !wasN) { // gap start
                
                pushbackGap = true;
                
                wasN = true;
                
            }else if (wasN && !isN(base)) { // internal gap end
                
                pushbackGap = true;
                
                wasN = false;
                
            }
            
            if(pos == (s.length()-1) && wasN && isN(base)) { // end of scaffold
                
                pushbackGap = true;
                pos++;
                
            }
            
            if (pushbackGap) {
                
                cout<<"g: "<<pos<<"\n";
                
                fastaGapBoundaries.push_back(pos);
                pushbackGap = false;
                
            }
            
            pos++;
            
        }
        
        setFastaGapBoundaries(fastaGapBoundaries);
        setFastaContigBoundaries(fastaGapBoundaries);
        
    }
    
    void setFastaHeader(string h) {
        fastaHeader = h;
    }
    
    void setFastaComment(string c) {
        fastaComment = c;
    }
    
    void setFastaSequence(string s) {
        fastaSequence = s;
    }
    
    void setFastaContigBoundaries(vector<unsigned int> fastaGapBoundaries) {
        
        vector<unsigned int> newFastaContigBoundaries;
        
        newFastaContigBoundaries.reserve(fastaGapBoundaries.size() + 2);
        
        if (fastaGapBoundaries.size() > 0) {
            
            newFastaContigBoundaries = fastaGapBoundaries;
            
            if (fastaGapBoundaries[0] != 0) {
                
                newFastaContigBoundaries.insert(newFastaContigBoundaries.begin(), 0);
                
            }else{
                
                newFastaContigBoundaries.erase(newFastaContigBoundaries.begin());
                
            }
            
            if (newFastaContigBoundaries[fastaGapBoundaries.size()] != fastaSequence.size()) {
                
                newFastaContigBoundaries.insert(newFastaContigBoundaries.end(), fastaSequence.size());
                
            }else{
                
                newFastaContigBoundaries.pop_back();
                
            }
            
        }
        
        else {
            
            newFastaContigBoundaries = {0, (unsigned int) fastaSequence.size()};
            
        }
        
        fastaContigBoundaries = newFastaContigBoundaries;
        
        for (unsigned int bound : newFastaContigBoundaries) {
            
            cout<<bound<<endl;
            
        }
        
    }
    
    void setFastaGapBoundaries(vector<unsigned int> g) {
        fastaGapBoundaries = g;
    }
    
    string getFastaHeader() {
        return fastaHeader;
    }
    
    string getFastaComment() {
        return fastaComment;
    }
    
    string getFastaSequence() {
        return fastaSequence;
    }
    
    unsigned int getFastaSeqLen() {
        return fastaSequence.size();
    }
    
    vector<unsigned int> getFastaContigs() {
        return fastaContigBoundaries;
    }
    
    vector<unsigned int> getFastaContigLens() {
        
        return bedIntervalSizes(fastaContigBoundaries);
    }
    
    vector<unsigned int> getFastaGaps() {
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
    
};

class FastaSequences {
    
private:
    vector<FastaSequence> newFasta = vector<FastaSequence>();
    
    vector<unsigned int> scaffLens;
    
    vector<unsigned int> contigLens;
    
    FastaSequence fastaSeq;
    
    unsigned long long int totScaffLen = 0;
    
    unsigned int totGapLen = 0, gapN = 0, scaffN50 = 0, scaffNG50 = 0, contigN50 = 0, contigNG50 = 0;
    
    double AverageScaffLen = 0;
    
    string h;
    char *c;
    
public:
    void appendFasta(string hg, string s) {
        
        h = string(strtok(strdup(hg.c_str())," ")); //process header line
        h.erase(0, 1);
        fastaSeq.setFastaHeader(h);
        
        c = strtok(NULL,""); //process comment line
        
        if (c != NULL) {
            
            fastaSeq.setFastaComment(string(c));
            
        }
        
        verbose(verbose_flag, "Header, comment, and fasta sequence read");
        
        verbose(verbose_flag, "Processing scaffold: "+h);
        
        fastaSeq.setFastaSequence(s);
        
        verbose(verbose_flag, "Fasta sequence set");
        
        fastaSeq.TraverseFastaSequence(s);
        
        verbose(verbose_flag, "Traversed fasta sequence");
        
        newFasta.push_back(fastaSeq);
        
        verbose(verbose_flag, "Fasta sequence added to fasta sequence vector");
        
        increaseTotScaffLen(fastaSeq.getFastaSeqLen());
        
        verbose(verbose_flag, "Increased total scaffold length");
        
        recordScaffLen(fastaSeq.getFastaSeqLen());
        
        verbose(verbose_flag, "Recorded length of fasta sequence");
        
        recordContigLens(fastaSeq.getFastaContigLens());
        
        verbose(verbose_flag, "Recorded length of contigs in fasta sequence");
        
        increaseTotGapLen(fastaSeq.getGapSum());
        
        verbose(verbose_flag, "Increased total gap length");
        
        increaseGapN(fastaSeq.getGapN());
        
        verbose(verbose_flag, "Increased total number of gaps");
        
        verbose(verbose_flag, "\n");
        
    }
    
    FastaSequence getFastaSequences(unsigned int idx) {
        
        FastaSequence fastaSequence = newFasta[idx];
        return fastaSequence;
        
    }
    
    unsigned int getScaffN() {
        
        return newFasta.size();
        
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
    
    void recordContigLens(vector <unsigned int> seqLens) {
        
        vector <unsigned int> newContigLens;
        
        newContigLens.reserve(contigLens.size() + seqLens.size());
        newContigLens.insert(newContigLens.end(), contigLens.begin(), contigLens.end());
        newContigLens.insert(newContigLens.end(), seqLens.begin(), seqLens.end());
        
        contigLens = newContigLens;
        
    }
    
    void computeScaffN50(unsigned int gSize) {
        
        sort(scaffLens.begin(), scaffLens.end(), greater<unsigned int>());
        
        unsigned long long int scaffSum = 0;
        
        for(unsigned int i = 0; i < getScaffN(); i++) {
            
            scaffSum += scaffLens[i];
            
            if (scaffSum >= getTotScaffLen() / 2 && scaffN50 < scaffLens[i]) {
                
                scaffN50 = scaffLens[i];
                
            }
            
            if (gSize > 0 && scaffSum >= gSize / 2 && scaffNG50 < scaffLens[i]) {
                
                scaffNG50 = scaffLens[i];
            }
            
            if (scaffN50 >= scaffLens[i] && scaffNG50 >= scaffLens[i]) {
                
                break;
                
            }
            
        }
        
    }
    
    void computeContigN50(unsigned int gSize) {
        
        sort(contigLens.begin(), contigLens.end(), greater<unsigned int>());
        
        unsigned long long int contigSum = 0;
        
        for(unsigned int i = 0; i < contigLens.size(); i++) {
            
            contigSum += contigLens[i];
            
            if (contigSum >= getTotScaffLen() / 2 && contigN50 < contigLens[i]) {
                
                contigN50 = contigLens[i];
                
            }
            
            if (gSize > 0 && contigSum >= gSize / 2 && contigNG50 < contigLens[i]) {
                
                contigNG50 = contigLens[i];
            }
            
            if (contigN50 >= contigLens[i] && contigNG50 >= contigLens[i]) {
                
                break;
                
            }
            
        }
        
    }
    
    unsigned int getScaffN50(unsigned long long int gSize) {
        
        computeScaffN50(gSize);
        
        return scaffN50;
        
    }
    
    unsigned int getScaffNG50() {
        
        return scaffNG50;
        
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
    
    unsigned int getLargestScaffold() {
        
        return scaffLens[0];
        
    }
    
    void setAverageScaffLen() {
        
        AverageScaffLen = (double) totScaffLen/scaffLens.size();
        
    }
    
    double getAverageScaffLen() {
        
        return AverageScaffLen;
        
    }
    
};

class FastaFile {
    
public:
    
    void ParseFasta(string newLine, FastaSequences &Fasta, string &fastaHeader, string &fastaSequence, unsigned int &idx) {
        
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
    
    
    FastaSequences Read(string iFileArg) {
        
        string newLine, fastaHeader, fastaSequence;
        unsigned int idx = 0;
        
        FastaSequences Fasta;
        
        ifstream stream(iFileArg);
        
        unsigned char buffer[2];
        stream.read((char*)(&buffer[0]), 2) ;
        
        
        stream.clear();
        stream.seekg(0, stream.beg);
        
        if (buffer[0] == 0x1f && (buffer[1] == 0x8b)) {
            
            stream.close();
            
            string fileData;
            if (!loadBinaryFile(iFileArg, fileData)) {
                printf("Error loading input file.");
            }
            
            string data;
            if (!gzipInflate(fileData, data)) {
                printf("Error decompressing input file.");
                
            }
            
            istringstream gzstream(data);
            
            while (getline (gzstream, newLine)) {
                
                if (gzstream) {
                    
                    ParseFasta(newLine, Fasta, fastaHeader, fastaSequence, idx);
                    
                }
                
                else {
                    
                    cout << "Gzip stream not successful.";
                    
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
