#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>

#include <unistd.h>
#include <getopt.h>

#include <vector>
#include <algorithm>

#include <chrono>

using namespace std;

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

static int tabular_flag;
string output(string output){
    
    if (tabular_flag) {
        
        output = output + "\t";
        
    }else{
        
        output = output + ": ";
        
    }
    
    return output;
}

class FastaSequence {
private:
    string fastaHeader;
    string fastaSequence;
    std::vector<int> fastaGaps;
    
public:
    
    void setFastaHeader(string h) {
        fastaHeader = h;
    }
    
    void setFastaSequence(string s) {
        fastaSequence = s;
    }
    
    void setFastaGaps(std::vector<int> g) {
        fastaGaps = g;
    }
    
    string getFastaSequence() {
        return fastaSequence;
    }
    
    
    string getFastaHeader() {
        return fastaHeader;
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
    std::vector<FastaSequence> newFasta = std::vector<FastaSequence>();
    
    FastaSequence fastaSeq;
    
public:
    void appendFasta(string h, string s) {
        h.erase(0, 1);
        fastaSeq.setFastaHeader(h);
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
    
};

class FastaFile {
    
public:
    FastaSequences Read(string iFileArg) {
        
        string newLine, fastaHeader, fastaSequence;
        int idx = 0;
        
        FastaSequences Fasta;
        
        ifstream thisFile(iFileArg);
        
        
        while (getline (thisFile, newLine)) {
            
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
        
        Fasta.appendFasta(fastaHeader,fastaSequence);
        
        thisFile.close();
        
        return Fasta;
        
    }
    
};