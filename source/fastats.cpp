//
//main.cpp
//xcode
//
//Created by Giulio Formenti on 12/17/21.
//

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
        
        int pos = 1;
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

int main(int argc, char **argv) {
    int c;
    int arg_counter;
    int digit_optind = 0;
    int pos_op = 1;
    int gSize = 0;
    
    string iFileArg;
    
    static int outSequence_flag;
    
    static int stats_flag;
    static int seqReport_flag;
    
    static int verbose_flag;
    static int cmd_flag;
    static int help_flag;
    
    if (argc == 1) {
        printf("in.fasta\n-h for additional help.\n");
        exit(0);
        
    }
    
    if (argc == 2) {
        stats_flag = 1; // default mode 'stats'
        
    }
    
    while (1) {
        
        int option_index = 0;
        
        static struct option long_options[] = {
            {"fasta", required_argument, 0, 'f'},
            
            {"out-sequence", no_argument, &outSequence_flag, 1},
            
            {"stats", no_argument, 0, 's'},
            {"seq-report", no_argument, &seqReport_flag, 1},
            {"tabular", no_argument, 0, 'v'},
            
            {"verbose", no_argument, 0, 'v'},
            {"cmd", no_argument, &cmd_flag, 1},
            {"help", no_argument, 0, 'h'},
            
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "-f:stvh",
                        long_options, &option_index);
        
        if (c == -1) {
            cout<<endl;
            break;
            
        }
        
        switch (c) {
            default:
                if (pos_op == 1) {iFileArg = optarg; pos_op++;}
                else if (pos_op == 2) {gSize = atoi(optarg); pos_op++;}
                else{printf("Error: too many positional arguments (%s).\n",optarg);exit(1);}
                break;
                
            case 0:
                break;
                
            case 'f':
                iFileArg = optarg;
                break;
                
            case 's':
                stats_flag = 1;
                break;
                
            case 't':
                tabular_flag = 1;
                break;
                
            case 'v':
                verbose_flag = 1;
                break;
                
            case 'h':
                printf("fastats in.fasta [genome size]\n");
                printf("Options:\n");
                printf("-f <file> fasta input. Also as first positional argument.\n");
                printf("-s report summary statistics.\n");
                printf("-t output in tabular format.\n");
                printf("-v verbose output.\n");
                exit(0);
        }
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    if (cmd_flag) {
        
        arg_counter = -1;
        while (arg_counter++ < argc-1) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
    FastaFile iFile;
    FastaSequences fastaSequences;
    
    fastaSequences = iFile.Read(iFileArg);
    int fastaSequencesN = fastaSequences.getScaffN();
    
    int counter = 0;
    FastaSequence fastaSequence;
    
    if (seqReport_flag) {
        
        while (counter < fastaSequencesN) {
            
            fastaSequence = fastaSequences.getFastaSequences(counter);
            
            cout<<"Seq "<<counter+1<<endl;
            cout<<"Header: "<<fastaSequence.getFastaHeader()<<endl;
            cout<<"Sequence length: "<<fastaSequence.getFastaSeqLen()<<endl;
            cout<<"Total gap length: "<<fastaSequence.gapSum()<<endl;
            cout<<"Number of Gaps: "<<fastaSequence.gapN()<<endl;
            
            if (outSequence_flag) {
                
                cout<<"Sequence: "<<fastaSequence.getFastaSequence()<<endl;
                
            }
            
            cout<<endl;
            counter++;
            
        }
        
        counter = 0;
        
    }
    
    if (stats_flag) {
        
        int totScaffLen = 0, scaffN50 = 0, scaffNG50 = 0, totGapLen = 0, gapN = 0;
        
        std::vector<int> scaffLens;
        
        while (counter < fastaSequencesN) {
            
            fastaSequence = fastaSequences.getFastaSequences(counter);
            scaffLens.push_back(fastaSequence.getFastaSeqLen());
            
            totScaffLen += scaffLens[counter];
            totGapLen += fastaSequence.gapSum();
            gapN += fastaSequence.gapN();
            
            counter++;
            
        }
        
        counter = 0;
        
        sort(scaffLens.begin(), scaffLens.end(), greater<int>());
        
        int scaffSum = 0, largestScaff = 0;
        
        for(int i = 0; i < fastaSequencesN; i++) {
            scaffSum += scaffLens[i];
            if (scaffSum >= totScaffLen / 2 && scaffN50 < scaffLens[i]) {
                
                scaffN50 = scaffLens[i];
                
            }
            
            if (gSize > 0 && scaffSum >= gSize / 2 && scaffNG50 < scaffLens[i]) {
                
                scaffNG50 = scaffLens[i];
            }
            
            if (scaffN50 >= scaffLens[i] && scaffNG50 >= scaffLens[i]) {
                
                break;
                
            }
            
        }
        
        cout<<output("N scaffold")<<fastaSequencesN<<endl;
        cout<<output("Total length")<<totScaffLen<<endl;
        cout<<output("Scaffold N50")<<scaffN50<<endl;
        
        if (gSize > 0) {
            
            cout<<output("Scaffold NG50")<<scaffNG50<<endl;
            
        }
        
        cout<<output("Largest scaffold")<<scaffLens[0]<<endl;
        cout<<output("Total gap length")<<totGapLen<<endl;
        cout<<output("Number of Gaps")<<gapN<<endl;
        
        counter = 0;
        
    }
    
    if (verbose_flag) {
        
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "\nElapsed time: " << elapsed.count() << " s\n";
        
        
    }
    
    exit(EXIT_SUCCESS);
    
}
