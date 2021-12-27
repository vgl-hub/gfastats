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

//functions
bool isN(char base){
    return (base == 'N' || base == 'n' || base == 'X');
}

int interval(std::vector<int> intervalVec, char op){
    
    int value = 0;
    
    if (!intervalVec.empty()){
        vector<int>::const_iterator end = intervalVec.cend();
        
        for (vector<int>::const_iterator it = intervalVec.cbegin(); it != end;) {
            
            switch (op) {
                    
                case 's': { // sum
                    
                    value += *(it+1) - *it;
                    break;
                }
                    
                case 'c': { // count
                    
                    value ++;
                    break;
                }
                    
            }
            
            it = it + 2;
        }
    }
    
    return value;
    
}

// string convertToString(char* a, char type)
// {
//
//     if (type){
//
//         s = "ss";
//
//     } else {
//
//         strig(a);
//
//     }
//
//     cout<<s;
//
//     return s;
// }

static int tabular_flag;
string output(string output){
    
    if (tabular_flag) {
        
        output = output + "\t";
        
    }else{
        
        output = output + ": ";
        
    }
    
    return output;
}

//zlib
bool gzipInflate( const std::string& compressedBytes, std::string& uncompressedBytes ) {
    if ( compressedBytes.size() == 0 ) {
        uncompressedBytes = compressedBytes ;
        return true ;
    }
    
    uncompressedBytes.clear() ;
    
    unsigned full_length = compressedBytes.size() ;
    unsigned half_length = compressedBytes.size() / 2;
    
    unsigned uncompLength = full_length ;
    char* uncomp = (char*) calloc( sizeof(char), uncompLength );
    
    z_stream strm;
    strm.next_in = (Bytef *) compressedBytes.c_str();
    strm.avail_in = compressedBytes.size() ;
    strm.total_out = 0;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    
    bool done = false ;
    
    if (inflateInit2(&strm, (16+MAX_WBITS)) != Z_OK) {
        free( uncomp );
        return false;
    }
    
    while (!done) {
        // If our output buffer is too small
        if (strm.total_out >= uncompLength ) {
            // Increase size of output buffer
            char* uncomp2 = (char*) calloc( sizeof(char), uncompLength + half_length );
            memcpy( uncomp2, uncomp, uncompLength );
            uncompLength += half_length ;
            free( uncomp );
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
        free( uncomp );
        return false;
    }
    
    for ( size_t i=0; i<strm.total_out; ++i ) {
        uncompressedBytes += uncomp[ i ];
    }
    free( uncomp );
    return true ;
}

/* Reads a file into memory. */
bool loadBinaryFile( const std::string& filename, std::string& contents ) {
    // Open the gzip file in binary mode
    FILE* f = fopen( filename.c_str(), "rb" );
    if ( f == NULL )
        return false ;
    
    // Clear existing bytes in output vector
    contents.clear();
    
    // Read all the bytes in the file
    int c = fgetc( f );
    while ( c != EOF ) {
        contents +=  (char) c ;
        c = fgetc( f );
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
    std::vector<int> fastaGaps;
    
public:
    
    void setFastaHeader(string h) {
        fastaHeader = h;
    }
    
    void setFastaComment(string c) {
        fastaComment = c;
    }
    
    void setFastaSequence(string s) {
        fastaSequence = s;
    }
    
    void setFastaGaps(std::vector<int> g) {
        fastaGaps = g;
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
    
    int getFastaSeqLen() {
        return fastaSequence.size();
    }
    
    void TraverseFastaSequence(string s) {
        
        unsigned int pos = 1;
        bool wasN = false;
        std::vector<int> fastaGaps;
        
        for (char base : s) {
            
            if (isN(base) && !wasN) { // gap start
                
                fastaGaps.push_back(pos);
                
                wasN = true;
                
            }else if(isN(base) && wasN && pos == s.length()) { // end of scaffold
                
                fastaGaps.push_back(pos+1);
                
            }
            
            if (!(isN(base)) && wasN) { // end of internal gap
                
                wasN = false;
                
                fastaGaps.push_back(pos);
            }
            
            pos++;
        }
        
        setFastaGaps(fastaGaps);
        
    }
    
    int gapSum() {
        
        return interval(fastaGaps,'s');
    }
    
    int gapN() {
        
        return interval(fastaGaps,'c');
    }
    
};

class FastaSequences {
    
private:
    vector<FastaSequence> newFasta = vector<FastaSequence>();
    
    vector<int> scaffLens;
    
    FastaSequence fastaSeq;
    
    long long int totScaffLen = 0;
    
    int totGapLen = 0, gapN = 0, scaffN50 = 0, scaffNG50 = 0;
    
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
        
        fastaSeq.setFastaSequence(s);
        fastaSeq.TraverseFastaSequence(s);
        
        newFasta.push_back(fastaSeq);
        
    }
    
    FastaSequence getFastaSequences(int idx) {
        
        FastaSequence fastaSequence = newFasta[idx];
        return fastaSequence;
        
    }
    
    int getScaffN() {
        
        return newFasta.size();
        
    }

    void increaseTotScaffLen(int ScaffLen) {

        totScaffLen += ScaffLen;
        
    }
    
    long long int getTotScaffLen() {

        return totScaffLen;
        
    }
    
    void increaseTotGapLen(int GapLen) {

        totGapLen += GapLen;
        
    }
    
    int getTotGapLen() {

        return totGapLen;
        
    }
    
    void increaseGapN(int GapN) {

        gapN += GapN;
        
    }
    
    int getTotGapN() {

        return gapN;
        
    }
    
    void recordScaffLen(int seqLen) {
    
        scaffLens.push_back(seqLen);
        
    }
    
    void computeScaffN50(int gSize, FastaSequences fastaSequences) {
        
        sort(scaffLens.begin(), scaffLens.end(), greater<int>());
        
        int scaffSum = 0;
        
        for(int i = 0; i < fastaSequences.getScaffN(); i++) {
            scaffSum += scaffLens[i];
            if (scaffSum >= fastaSequences.getTotScaffLen() / 2 && scaffN50 < scaffLens[i]) {
                
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
    
    int getScaffN50() {

        return scaffN50;
        
    }
    
    int getScaffNG50() {

        return scaffNG50;
        
    }
    
    int getLargestScaffold() {

        return scaffLens[0];
        
    }
    
};

class FastaFile {
    
public:
    
    void ParseFasta(string newLine, FastaSequences &Fasta, string &fastaHeader, string &fastaSequence, int &idx) {
        
        switch (newLine[0]) {
                
            case '>': {
                
                if (idx> 0 ) {
                    
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
        int idx = 0;
        
        FastaSequences Fasta;
        
        ifstream stream(iFileArg);
        
        unsigned char buffer[2];
        stream.read( (char*)(&buffer[0]), 2 ) ;
        
        
        stream.clear();
        stream.seekg(0, stream.beg);
        
        if (buffer[0] == 0x1f && (buffer[1] == 0x8b)) {
            
            stream.close();
            
            string fileData;
            if (!loadBinaryFile(iFileArg, fileData)) {
                printf( "Error loading input file." );
            }
            
            string data;
            if (!gzipInflate(fileData, data)) {
                printf( "Error decompressing input file." );
                
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
