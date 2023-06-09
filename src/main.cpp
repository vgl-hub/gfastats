//
//fastats.cpp
//
//Created by Giulio Formenti on 12/17/21.
//

#include "main.h"

std::string version = "1.3.6";

//global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

short int tabular_flag;
int verbose_flag;
int seqReport_flag;
int outSequence_flag;
int nstarReport_flag;
int outSize_flag;
int outCoord_flag;
int outFile_flag;
int outBubbles_flag;
int stats_flag;
int cmd_flag;
int rmGaps_flag;
int discoverPaths_flag;
int extractContigs_flag;
int hc_flag;
int hc_cutoff;
int terminalOvlLen = 0;
int maxThreads = 0;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    unsigned long long int gSize = 0; // expected genome size, with 0 NG/LG* statistics are not computed
    int splitLength = 0; // line length for fasta output
    
    char bedOutType = 'a'; // default output type with bed flag (agp)
    
    bool isPipe = false; // to check if input is from pipe
    
    UserInputGfastats userInput; // initialize input object
    
    if (argc == 1) { // gfastats with no arguments
            
        printf("gfastats input.[fasta|fastq|gfa][.gz] [expected genome size] [header[:start-end]]\n-h for additional help.\n");
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        
        {"threads", required_argument, 0, 'j'},
        
        {"agp-to-path", required_argument, 0, 'a'}, // agp to path conversion
        {"swiss-army-knife", required_argument, 0, 'k'}, // the swiss army knife
        {"remove-terminal-gaps", no_argument, &rmGaps_flag, 1}, // this remove all gap edges at the end of sequences
        {"homopolymer-compress", required_argument, &hc_flag, 1},
        {"discover-paths", no_argument, &discoverPaths_flag, 1},
        {"discover-terminal-overlaps", optional_argument, 0, 0},
        {"sort", required_argument, 0, 0},
        {"extract-contigs", no_argument, &extractContigs_flag, 1},
        
        {"include-bed", required_argument, 0, 'i'},
        {"exclude-bed", required_argument, 0, 'e'},
 
        {"out-format", required_argument, 0, 'o'},
        {"line-length", required_argument, 0, 0},
        {"no-sequence", no_argument, &userInput.noSequence, 1},
        {"out-sequence", no_argument, &outSequence_flag, 1},
        {"out-size", required_argument, 0, 's'},
        {"out-coord", required_argument, 0, 'b'},
        {"out-bubbles", no_argument, &outBubbles_flag, 1},
        
        {"stats", no_argument, &stats_flag, 1},
        {"seq-report", no_argument, &seqReport_flag, 1},
        {"nstar-report", no_argument, &nstarReport_flag, 1},
        {"tabular", no_argument, 0, 't'},
        {"locale", required_argument, 0, 0},
        
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    
    while (1) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:a:b:e:f:i:j:k:o:s:tvh",
                        long_options, &option_index);

        if (optind < argc && !isPipe) { // if pipe wasn't assigned already
            
            isPipe = isDash(argv[optind]) ? true : false; // check if the argument to the option is a '-' and set it as pipe input
            
        }
        
        if (optarg != nullptr && !isPipe) { // case where pipe input is given as positional argument (input sequence file)
        
            isPipe = isDash(optarg) ? true : false;
            
        }

        if (c == -1) { // exit the loop if run out of options
            break;
            
        }
        
        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    case 'b':
                        bedOutType = 'a'; // default bed output is agp is -b option is given without argument
                        outCoord_flag = 1;
                        break;
                        
                    case 's':
                        bedOutType = 's'; // default size output is scaffold is -s option is given without argument
                        outSize_flag = 1;
                        break;
                        
                    case 'o':
                        userInput.outFiles.push_back("fasta"); // default output is fasta if -o option is given without argument
                        outFile_flag = 1;
                        break;
                        
                    default:
                        fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                        return EXIT_FAILURE;
                }
                break;
            default: // handle positional arguments
                if (pos_op == 1) { // only one positional argument given
                    
                    if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                    
                        userInput.pipeType = 'f'; // pipe input is a sequence
                    
                    }else{ // input is a regular file
                        
                        ifFileExists(optarg);
                        userInput.inSequence = optarg;
                        
                    }
                    
                    pos_op++;
                    
                }else if (pos_op == 2 || pos_op == 3) { // if >1 positional argument, check what additional positional arguments are present
                    
                    if (isInt(optarg)) { // if the positional argument is a number, it is likely the expected genome size
                        
                        gSize = atoll(optarg); pos_op++;
                        
                    }else{ // else it is an include argument
                        
                        std::string header = optarg, cBegin, cEnd; // the header for coordinates provided as positional argument
                        
                        reverse(header.begin(), header.end()); // we work our way from the end
                        
                        cBegin = header.substr(header.find('-') + 1, header.find(':') - header.find('-') - 1);
                        cEnd = header.substr(0, header.find('-'));
                        
                        if(isNumber(cEnd) &&
                           isNumber(cBegin)) { // prevent headers with : - characters to break the extraction
                            
                            header = header.substr(header.find(':') + 1, header.size());
                            reverse(header.begin(), header.end());
                            
                            reverse(cBegin.begin(), cBegin.end());
                            reverse(cEnd.begin(), cEnd.end());
                            
                        }else{
                            
                            reverse(header.begin(), header.end());
                            cBegin = "0";
                            cEnd = "0";
                            
                        }
                        
                        userInput.bedIncludeList.pushCoordinates(header, stoi(cBegin), stoi(cEnd)); pos_op++;
                        
                    }
                    
                }
                else{printf("Error: too many positional arguments (%s).\n",optarg);exit(1);}
                
                break;
                
            case 0: // case for long options without short options
                
                if (strcmp(long_options[option_index].name,"discover-terminal-overlaps") == 0) {
                    
                    if (optarg == NULL && optind < argc
                        && argv[optind][0] != '-')
                    {
                        optarg = argv[optind++];
                    }
                    
                    if (optarg != NULL) {
                        
                        terminalOvlLen = atoi(optarg);
                        
                    }else {
                        
                        terminalOvlLen = 1000;
                        
                    }
                    
                }

                if (strcmp(long_options[option_index].name,"line-length") == 0)
                    splitLength = atoi(optarg);
                
                if (strcmp(long_options[option_index].name,"sort") == 0) {

                    std::vector<std::string> options {"none", "ascending", "descending", "largest", "smallest"};
                            
                    if (std::find(options.begin(), options.end(), optarg) != options.end() || ifFileExists(optarg)){
                        
                        userInput.sortType = optarg;
                        
                    }else{printf("Error: unrecognized sorting option (%s).\n", optarg); exit(1);}
                    
                }

                if(strcmp(long_options[option_index].name,"homopolymer-compress") == 0) {
                    hc_cutoff = atoi(optarg);
                    stats_flag = true;
                }
                
                if (strcmp(long_options[option_index].name,"locale") == 0) {
                    
                    setlocale(LC_ALL, optarg);
                    std::cout.imbue(std::locale(optarg));
                    std::locale::global(std::locale(optarg));
                    stats_flag = true;
                    
                }
                
                break;
                
            case 'a': // agp to paths
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'a'; // pipe input is agp
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.inAgp = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
            
            case 'b': // output bed type (agp, contig, gaps)
                bedOutType = *optarg;
                outCoord_flag = 1;
                break;
                
            case 'e': // bed exclude
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'e'; // pipe input is an exclude bed
                
                }else{ // input is a regular file
                
                    ifFileExists(optarg);
                    userInput.inBedExclude = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
                
            case 'f': // input sequence
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'f'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.inSequence = optarg;
                    stats_flag = true;
                    
                }
                    
                break;
                
            case 'i': // bed include
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'i'; // pipe input is an include bed
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.inBedInclude = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
                
            case 'j': // max threads
                maxThreads = atoi(optarg);
                stats_flag = 1;
                break;
                
            case 'k': // the swiss army knife
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'k'; // pipe input is a set of instructions for the swiss army knife
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.inSak = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
                
            case 'o': // handle output (file or stdout)
                
                outFile_flag = 1;

                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'r'; // pipe input is a sequence
                    
                }else{ // outputs are regular files
                    
                    optind--;
                    
                    std::string file;
                    uint8_t i = 0;

                    for( ;optind < argc && !isInt(argv[optind]); optind++) {
                        
                        if (i > 0 && *argv[optind] == '-')
                            break;
                        
                        file = argv[optind];
                        
                        if (file.find("-o") != std::string::npos)
                            file.erase(0, 2); // handle file name attached to option
                        
                        userInput.outFiles.push_back(file);
                        
                        ++i;
                        
                    }
                    
                    stats_flag = true;
                    
                }
                
                break;
                
            case 's': // output size of features
                bedOutType = *optarg;
                outSize_flag = 1;
                break;
                
            case 't': // tabular output
                tabular_flag = 1;
                break;
                
            case 'v': // software version
                printf("gfastats v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                exit(0);
                
            case 'h': // help
                printf("gfastats input.[fasta|fastq|gfa][.gz] [expected genome size] [header[:start-end]]\n");
                printf("genome size: estimated genome size for NG* statistics (optional).\n");
                printf("header: target specific sequence by header, optionally with coordinates (optional).\n");
                printf("\nOptions:\n");
                printf("\t-a --agp-to-path <file> converts input agp to path and replaces existing paths.\n");
                printf("\t-b --out-coord a|s|c|g generates bed coordinates of given feature (agp|scaffolds|contigs|gaps default:agp).\n");
                printf("\t-e --exclude-bed <file> opposite of --include-bed. They can be combined (no coordinates).\n");
                printf("\t-f --input-sequence <file> input file (fasta, fastq, gfa [.gz]). Also as first positional argument.\n");
                printf("\t-h --help print help and exit.\n");
                printf("\t-i --include-bed <file> generates output on a subset list of headers or coordinates in 0-based bed format.\n");
                printf("\t-k --swiss-army-knife <file> set of instructions provided as an ordered list.\n");
                printf("\t-j --threads <n> numbers of threads (default: max).\n");
                printf("\t-o --out-format fasta|fastq|gfa[.gz] outputs selected sequences. If more than the extension is provided the output is written to the specified file (e.g. out.fasta.gz). Multiple file outputs can be given at once.\n");
                printf("\t-s --out-size s|c|g  generates size list of given feature (scaffolds|contigs|gaps default:scaffolds).\n");
                printf("\t-t --tabular output in tabular format.\n");
                printf("\t-v --version software version.\n\n");
                printf("\t--cmd print $0 to stdout.\n");
                printf("\t--discover-paths prototype to induce paths from input.\n");
                printf("\t--discover-terminal-overlaps <n> append perfect terminal overlaps of minimum length n (default: 1000).\n");
                printf("\t--homopolymer-compress <n> compress all the homopolymers longer than n in the input.\n");
                printf("\t--line-length <n> specifies line length in when output format is fasta. Default has no line breaks.\n");
                printf("\t--nstar-report generates full N* and L* statistics.\n");
                printf("\t--out-sequence reports also the actual sequence (in combination with --seq-report).\n");
                printf("\t--out-bubbles outputs a potential list of bubbles in the graph.\n");
                printf("\t--seq-report report statistics for each sequence.\n");
                printf("\t--sort ascending|descending|largest|smallest|file sort sequences according to input. Ascending/descending used the sequence/path header.\n");
                printf("\t--stats report summary statistics (default).\n");
                printf("\t--verbose verbose output.\n");
                printf("\t--locale set a different locale, for instance to use , for thousand separators use en_US.UTF-8.\n");
                printf("\nAll input files can be piped from stdin using '-'.\n");
                exit(0);
        }
        
        if    (argc == 2 || // handle various cases in which the output should include summary stats
              (argc == 3 && pos_op == 2) ||
              (argc == 4 && pos_op == 3) ||
               nstarReport_flag ||
               discoverPaths_flag) {
            
            stats_flag = 1; // default mode 'stats'
            
        }
        
    }
    
    lg.verbose("Input variables assigned");
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
    Input in;
    
    in.load(userInput); // load user input
    
    lg.verbose("Loaded user input");
    
    InSequences inSequences; // initialize sequence collection object
    
    lg.verbose("Sequence object generated");
    
    in.read(inSequences); // read input content to inSequences container

    if(verbose_flag) {std::cerr<<"\n\n";};

    std::vector<Log> logs = inSequences.getLogs();

    //consolidate log
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    lg.verbose("Finished reading input files");
    if(verbose_flag) {std::cerr<<"\n";};
    
    Report report;
    
    if (seqReport_flag || outSequence_flag) { // report results for each sequence
        
        stats_flag = false;
        
        report.seqReport(inSequences, outSequence_flag);
        
    }
    
    if (outFile_flag) { // output sequences to file or stdout
        
        stats_flag = false;
        
        for (std::string file : userInput.outFiles)
            report.outFile(inSequences, file, userInput, splitLength);
        
    }
    
    if (outCoord_flag || outSize_flag) { // output coordinates
        
        stats_flag = false;
        
        report.outCoord(inSequences, bedOutType, outSize_flag);
        
    }
    
    if (stats_flag) { // output summary statistics
        
        report.reportStats(inSequences, gSize);
        
    }
    
    if (nstarReport_flag) { // output full N/L* statistics
        
        report.nstarReport(inSequences, gSize);
        
    }
    
    lg.verbose("Generated output");
    
    exit(EXIT_SUCCESS);
    
}
