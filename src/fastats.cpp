//
//fastats.cpp
//xcode
//
//Created by Giulio Formenti on 12/17/21.
//

//global
static int verbose_flag;
static int seqReport_flag;

#include <fastats.h>

int main(int argc, char **argv) {
    int c;
    int arg_counter;
    int pos_op = 1;
    int gSize = 0;
    
    string iFileArg;
    
    static int outSequence_flag;
    
    static int stats_flag;

    static int cmd_flag;
    
    if (argc == 1) {
        printf("in.fasta\n-h for additional help.\n");
        exit(0);
        
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
                printf("-f --fasta <file> fasta input. Also as first positional argument.\n");
                printf("-s --stats report summary statistics.\n");
                printf("-t --tabular output in tabular format.\n");
                printf("-v --verbose verbose output.\n");
                printf("-h --help print help and exit.\n");
                printf("--seq-report report statistics for each sequence.\n");
                printf("--out-sequence reports also the actual sequence (in combination with --seq-report).\n");                
                printf("--cmd print $0 to stdout.\n");
                exit(0);
        }
        
        if (argc == 2 || (argc == 3 && pos_op ==2)) {
        
        	stats_flag = 1; // default mode 'stats'
        
    	}
        
    }
    
    if (cmd_flag) {
        
        arg_counter = -1;
        while (arg_counter++ < argc-1) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
    FastaFile iFile;
    
    verbose(verbose_flag, "File object generated");
    
    FastaSequences fastaSequences;

    verbose(verbose_flag, "Fasta sequence object generated");
    
    fastaSequences = iFile.Read(iFileArg);
    
    verbose(verbose_flag, "Finished reading sequences from file to fasta sequence object");
    
    int counter = 0;
    FastaSequence fastaSequence;
    
    if (seqReport_flag) {
        
        while (counter < fastaSequences.getScaffN()) {
            
            fastaSequence = fastaSequences.getFastaSequences(counter);
            
            cout<<"Seq "<<counter+1<<endl;
            cout<<"Header: "<<fastaSequence.getFastaHeader()<<endl;
            cout<<"Comment: "<<fastaSequence.getFastaComment()<<endl;
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
        
        fastaSequences.computeScaffN50(gSize, fastaSequences);
        
        verbose(verbose_flag, "Computed scaffN50");
        
        cout<<output("N scaffold")<<fastaSequences.getScaffN()<<endl;
        cout<<output("Total length")<<fastaSequences.getTotScaffLen()<<endl;
        cout<<output("Scaffold N50")<<fastaSequences.getScaffN50()<<endl;
        
        if (gSize > 0) {
            
            cout<<output("Scaffold NG50")<<fastaSequences.getScaffNG50()<<endl;
            
        }
        
        cout<<output("Largest scaffold")<<fastaSequences.getLargestScaffold()<<endl;
        cout<<output("Total gap length")<<fastaSequences.getTotGapLen()<<endl;
        cout<<output("Number of Gaps")<<fastaSequences.getTotGapN()<<endl;
        
        counter = 0;
        
    }
    
    verbose(verbose_flag, "Generated output");
    
    exit(EXIT_SUCCESS);
    
}
