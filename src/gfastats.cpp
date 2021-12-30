//
//fastats.cpp
//
//Created by Giulio Formenti on 12/17/21.
//

//global
static int verbose_flag;
static int seqReport_flag;

#include <gfastats.h>

int main(int argc, char **argv) {
    
    short int c;
    short unsigned int arg_counter;
    short unsigned int pos_op = 1;
    unsigned long long int gSize = 0;
    
    std::string iFileArg;
    
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
            {"tabular", no_argument, 0, 't'},
            
            {"verbose", no_argument, 0, 'v'},
            {"cmd", no_argument, &cmd_flag, 1},
            {"help", no_argument, 0, 'h'},
            
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "-f:stvh",
                        long_options, &option_index);
        
        if (c == -1) {
            std::cout<<std::endl;
            break;
            
        }
        
        switch (c) {
            default:
                if (pos_op == 1) {iFileArg = optarg; pos_op++;}
                else if (pos_op == 2) {gSize = atoi(optarg); pos_op++;}
                else{printf("Error: too many positional arguments (%s).\n",optarg);exit(1);}
                
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
    
    unsigned int counter = 0;
    FastaSequence fastaSequence;
    
    if (seqReport_flag) {
        
        while (counter < fastaSequences.getScaffN()) {
            
            fastaSequence = fastaSequences.getFastaSequences(counter);
            
            std::cout<<output("Seq:")<<counter+1<<std::endl;
            std::cout<<output("Header:")<<fastaSequence.getFastaHeader()<<std::endl;
            std::cout<<output("Comment:")<<fastaSequence.getFastaComment()<<std::endl;
            std::cout<<output("Sequence length:")<<fastaSequence.getFastaSeqLen()<<std::endl;
            std::cout<<output("Total gap length:")<<fastaSequence.getGapSum()<<std::endl;
            std::cout<<output("Number of Gaps:")<<fastaSequence.getGapN()<<std::endl;

            printf("%s%u, %u, %u, %u\n",output("Base composition (ACGT):").c_str(), fastaSequence.getA(),
                fastaSequence.getC(),
                fastaSequence.getG(),
                fastaSequence.getT());
            printf("%s%.2f\n",output("GC content %:").c_str(), fastaSequence.computeGCcontent());

            
            if (outSequence_flag) {
                
                std::cout<<output("Sequence:")<<fastaSequence.getFastaSequence()<<std::endl;
                
            }
            
            std::cout<<std::endl;
            counter++;
            
        }
        
        counter = 0;
        
        std::cout<<output("+++Summary+++")<<std::endl;
        
    }
    
    if (stats_flag) {
        
        verbose(verbose_flag, "Computed scaffN50");
        
        std::cout<<output("N scaffolds:")<<fastaSequences.getScaffN()<<std::endl;
        std::cout<<output("Total length:")<<fastaSequences.getTotScaffLen()<<std::endl;
        printf("%s%.2f\n",output("Average scaffold length:").c_str(), fastaSequences.computeAverageScaffLen());
        std::cout<<output("Scaffold N50:")<<fastaSequences.getScaffN50(gSize)<<std::endl;
        std::cout<<output("Scaffold L50:")<<fastaSequences.getScaffL50()<<std::endl;
        
        if (gSize > 0) {
            
            std::cout<<output("Scaffold NG50:")<<fastaSequences.getScaffNG50()<<std::endl;
            std::cout<<output("Scaffold LG50:")<<fastaSequences.getScaffLG50()<<std::endl;
            
        }
        std::cout<<output("N contigs:")<<fastaSequences.getContigN()<<std::endl;
        std::cout<<output("Contig N50:")<<fastaSequences.getContigN50(gSize)<<std::endl;
        std::cout<<output("Contig L50:")<<fastaSequences.getContigL50()<<std::endl;
        
        if (gSize > 0) {
            
            std::cout<<output("Contig NG50:")<<fastaSequences.getContigNG50()<<std::endl;
            std::cout<<output("Contig LG50:")<<fastaSequences.getContigLG50()<<std::endl;
            
        }
        std::cout<<output("Largest scaffold:")<<fastaSequences.getLargestScaffold()<<std::endl;
        std::cout<<output("Total gap length:")<<fastaSequences.getTotGapLen()<<std::endl;
        std::cout<<output("Number of Gaps:")<<fastaSequences.getTotGapN()<<std::endl;

        printf("%s%lu, %lu, %lu, %lu\n",output("Base composition (ACGT):").c_str(), fastaSequences.getTotA(),
            fastaSequences.getTotC(),
            fastaSequences.getTotG(),
            fastaSequences.getTotT());
        printf("%s%.2f\n",output("GC content %:").c_str(), fastaSequences.computeGCcontent());
        
        counter = 0;
        
    }
    
    verbose(verbose_flag, "Generated output");
    
    exit(EXIT_SUCCESS);
    
}
