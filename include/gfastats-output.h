//
//  gfastats-classes.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef gfastatsoutput_h
#define gfastatsoutput_h

//classes
class Report {
private:
    unsigned int counter = 0;
    
public:
    bool seqReport(InSequences &inSequences, InSegment &inSegment, int &outSequence_flag) { // method to output the summary statistics for each sequence
        
        counter = 0;

        while (counter < inSequences.getSegmentN()) {

            inSegment = inSequences.getInSegment(counter);

            std::cout<<output("Seq")<<counter+1<<"\n";
            std::cout<<output("Header")<<inSegment.getSeqHeader()<<"\n";
            std::cout<<output("Comment")<<inSegment.getSeqComment()<<"\n";
            std::cout<<output("Total segment length")<<inSegment.getSegmentLen()<<"\n";

            printf("%s%u, %u, %u, %u\n",output("Base composition (ACGT)").c_str(), inSegment.getA(),
                   inSegment.getC(),
                   inSegment.getG(),
                   inSegment.getT());
            printf("%s%.2f\n",output("GC content %").c_str(), inSegment.computeGCcontent());
            std::cout<<output("# soft-masked bases")<<inSegment.getLowerCount()<<"\n";


            if (outSequence_flag) {

                std::cout<<output("Sequence")<<inSegment.getInSequence()<<"\n";
                std::cout<<output("Quality")<<inSegment.getInSequenceQuality()<<"\n";

            }

            std::cout<<"\n";
            counter++;

        }

        counter = 0;
        
        return true;
//          for the scaff report
//        std::cout<<output("# contig")<<inSegment.getContigN()<<"\n";
//        std::cout<<output("Total sequence length")<<inSegment.getSegmentLength()<<"\n";
//        std::cout<<output("Total gap length")<<inSegment.getGapSum()<<"\n";
//        std::cout<<output("# gaps")<<inSegment.getGapN()<<"\n";
        
    }
    
    bool outFile(InSequences &inSequences, InSegment &inSegment, int splitLength, std::string &outSeq) { // method to output new sequence opposed to sequence report
        
        // unordered map to handle out correspondence in following switch statement
        const static std::unordered_map<std::string,int> string_to_case{
            {"fasta",1},
            {"fa",1},
            {"fasta.gz",1},
            {"fa.gz",1},
            {"fastq",2},
            {"fq",2},
            {"fastq.gz",2},
            {"fq.gz",2},
            {"gfa",3},
            {"gfa.gz",3}
        };
        
        // variables to handle output type
        bool gzip = false;
        bool outFile = false;
        
        // variable to handle output path and extension
        std::string path = rmFileExt(outSeq);
        std::string ext = getFileExt(outSeq);
        
        // depending on use input get output format
        if(ext == "gz") {
            
            ext = getFileExt(path) + ".gz";
            path = rmFileExt(path);
            gzip = true;
            
        }
        
        // if the input is not in the unordered map, it means we need to write a new file with the path provided by the user otherwise the output is in the format specified by the user
        if (string_to_case.find(path) == string_to_case.end()) {
            
            outFile = true;
            
        }else{
            
            ext = outSeq;
            
        }
        
        // here we create a smart pointer to handle any kind of output stream
        std::unique_ptr<std::ostream> stream;
        
        // this stream outputs to file
        std::ofstream ofs(outSeq);

        if (gzip && outFile) { // if the requested output is gzip compressed and should be outputted to a file
            
            // this stream outputs gzip compressed to file
            zstream::ogzstream zfout(ofs);
            
            stream = make_unique<std::ostream>(zfout.rdbuf()); // then we use the stream for gzip compressed file outputs
            
        }else if (!gzip && outFile){ // else if no compression is requested
            
            stream = make_unique<std::ostream>(ofs.rdbuf());  // we use the stream regular file outputs
            
        }else{ // else the output is not written to a file
            
            // we close and delete the file
            ofs.close();
            remove(outSeq.c_str());
            
            if (gzip) { // if the output to stdout needs to be compressed we use the appropriate stream
            
                // this stream outputs gzip compressed to stdout
                zstream::ogzstream zout(std::cout);
                
                stream = make_unique<std::ostream>(zout.rdbuf());
            
            }else{ // else we use a regular cout stream
                
                stream = make_unique<std::ostream>(std::cout.rdbuf());
                
            }
        }
        
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) { // this switch allows us to generate the output according to the input request and the unordered map. If the requested output format is not in the map we fall back to the undefined case (0)
                
            case 1: { // fasta[.gz]
                
                std::string seqHeader;
                
                std::string inSeq; // the new sequence being built recursively
                
                // generate adjacency list representation of a graph
                inSequences.buildGraph(inSequences.getGFAGaps());
                
                for (unsigned int i = 0; i != inSequences.getAdjListFW().size(); ++i) { // loop through all nodes
                    
                    InSegment inSegment; // a new inSequence object, the result of concatenating by gaps
                    
                    if (!inSequences.getVisited(i)) { // check if the node was already visited
                        
                        verbose(verbose_flag, "Graph DFS");
                        inSequences.dfsSeq(i, inSeq); // if not, visit all connected components recursively
                        
                        seqHeader = inSequences.getInSegment(i).getSeqHeader();
                        
                        seqHeader.pop_back();
                        seqHeader.pop_back();
                        
                        *stream<<">"<<seqHeader<<" "<<inSequences.getInSegment(i).getSeqComment()<<"\n";
                        
                        if (splitLength != 0) {
                            
                            textWrap(inSeq, *stream, splitLength); // wrapping text at user-specified line length
                            
                        }
                        
                        *stream<<inSeq<<"\n";
                        inSeq = "";
                        
                    }
                    
                }
                
                break;
                
            }
                
            case 2: { // fastq[.gz]
                
                std::string inSeq, inSeqQual; // the new sequence being built recursively and its quality
                
                // generate adjacency list representation of a graph
                inSequences.buildGraph(inSequences.getGFAGaps());
                
                for (unsigned int i = 0; i != inSequences.getAdjListFW().size(); ++i) { // loop through all nodes
                    
                    InSegment inSegment; // a new inSequence object, the result of concatenating by gaps
                    
                    if (!inSequences.getVisited(i)) { // check if the node was already visited
                        
                        verbose(verbose_flag, "Graph DFS");
                        inSequences.dfsSeq(i, inSeq, &inSeqQual); // if not, visit all connected components recursively
                        
                        *stream<<"@"<<inSequences.getInSegment(i).getSeqHeader()<<" "<<inSequences.getInSegment(i).getSeqComment()<<"\n"<<inSeq<<"\n+\n"<<inSeqQual<<"\n";
                        
                        inSeq = "";
                        inSeqQual = "";
                        
                    }
                    
                }
                
                break;
                
            }
                
            case 3: { // gfa[.gz]
                
                std::string seqHeader;
                
                *stream<<"H\tVN:Z:2.0\n";
                
                for (InSegment inSegment : inSequences.getInSegments()) {
                    
                    seqHeader = inSegment.getSeqHeader();
                    
                    *stream <<"S\t" // line type
                            <<seqHeader<<"\t" // header
                            <<inSegment.getSegmentLen()<<"\t" // seq length
                            <<inSegment.getInSequence(); // sequence
                    
                    if (inSegment.getSeqComment() != "") {
                        
                        *stream <<"\tC:"<<inSegment.getSeqComment(); // optional comment
                        
                    }
                    
                    if (inSegment.getInSequenceQuality() != "") {
                        
                        *stream <<"\tQ:"<<inSegment.getInSequenceQuality(); // optional comment
                        
                    }
                    
                    *stream<<"\n";
                    
                }
                
                for (InGap inGap : inSequences.getGFAGaps()) {
                    
                    *stream <<"G\t" // line type
                            <<inGap.getgId()<<"\t" // id
                            <<inGap.getsId1()<<inGap.getsId1Or()<<"\t" // sid1:ref
                            <<inGap.getsId2()<<inGap.getsId2Or()<<"\t" // sid2:ref
                            <<inGap.getDist()<<"\n"; // size
                    
                }
                
                break;
                
            }
                
            case 0: { // undefined case
                
                std::cout<<"Unrecognized output: "<<outSeq;
                
                break;
                
            }
                
        }
        
        if(outFile) { // if we wrote to file, we close it
            
            ofs.close();
            
        }
        
        return true;
        
    }
    
    bool outSize(InSequences &inSequences, InSegment &inSegment, char &sizeOutType) { // method to output only the size of the sequences
        
        std::string seqHeader;
        std::vector<unsigned int> seqBoundaries;
        
        switch (sizeOutType) {
 
            default:
            case 's': { // scaffolds
                
                std::string seqHeader, seqComment, inSeq; // header, comment and the new sequence being built recursively
                
                // generate adjacency list representation of a graph
                inSequences.buildGraph(inSequences.getGFAGaps());
                
                for (unsigned int i = 0; i != inSequences.getAdjListFW().size(); ++i) { // loop through all nodes
                    
                    InSegment inSegment; // a new inSequence object, the result of concatenating by gaps
                
                    seqHeader = inSequences.getInSegment(i).getSeqHeader();
                    seqComment = inSequences.getInSegment(i).getSeqComment();
                    
                    if (!inSequences.getVisited(i)) { // check if the node was already visited
                        
                        verbose(verbose_flag, "Graph DFS");
                        inSequences.dfsSeq(i, inSeq); // if not, visit all connected components recursively
                        
                        seqHeader.pop_back();
                        seqHeader.pop_back();
                        
                        std::cout<<seqHeader<<"\t"<<inSeq.size()<<"\n";
                        
                        inSeq = "";
                        
                    }
                    
                }
                
                break;
            }
                
            case 'c': { // contigs
                
                for (InSegment inSegment : inSequences.getInSegments()) {
                    
                    seqHeader = inSegment.getSeqHeader(); // remove unique internal identifier
                    seqHeader.pop_back();
                    seqHeader.pop_back();
                    
                    std::cout<<seqHeader<<"\t"<<inSegment.getInSequence().size()<<"\n";

                }
                
                break;
                
            }
                
            case 'g': { // gaps
                
                for (InGap inGap : inSequences.getInGaps()) {
                    
                    seqHeader = inGap.getgId(); // remove unique internal identifier
                    seqHeader.pop_back();
                    seqHeader.pop_back();

                    std::cout<<seqHeader<<"\t"<<inGap.getDist()<<"\n";

                }
                
                break;
                
            }
                
        }
        
        return true;
        
    }
    
    bool outCoord(InSequences &inSequences, InSegment &inSegment, char bedOutType) { // method to output the coordinates of each feature
        
        counter = 0;
        
        std::string seqHeader;
        std::vector<unsigned int> seqBoundaries;
        
        switch (bedOutType) {
                
            case 'c': { // contigs
                
                unsigned int pos = 0, gapLen = 0;
                bool wasN = false;
                        
                std::string seqHeader, seqComment, inSeq; // header, comment and the new sequence being built recursively
                
                // generate adjacency list representation of a graph
                inSequences.buildGraph(inSequences.getGFAGaps());
                
                for (unsigned int i = 0; i != inSequences.getAdjListFW().size(); ++i) { // loop through all nodes
                    
                    InSegment inSegment; // a new inSequence object, the result of concatenating by gaps
                
                    seqHeader = inSequences.getInSegment(i).getSeqHeader();
                    seqComment = inSequences.getInSegment(i).getSeqComment();
                    
                    if (!inSequences.getVisited(i)) { // check if the node was already visited
                        
                        verbose(verbose_flag, "Graph DFS");
                        inSequences.dfsSeq(i, inSeq); // if not, visit all connected components recursively
                        
                        seqHeader.pop_back();
                        seqHeader.pop_back();
                        
                        for (char &base : inSeq) {

                            switch (base) {

                                case 'N': {
                                    
                                    gapLen++;
                                    
                                    if (!wasN) { // segment end
                                    
                                        std::cout<<"\t"<<pos<<"\n";
                                        
                                    }
                                    
                                    wasN = true;
                                    
                                    break;
                                    
                                }
                                default: {
                                    
                                    if (pos == 0) { // sequence start

                                        std::cout<<seqHeader<<"\t"<<0;

                                    }
                                    
                                    if (wasN) { // segment start
                                        
                                        std::cout<<seqHeader<<"\t"<<pos;
                                        
                                        
                                    }
                                    
                                    if (pos == inSeq.size()-1) { // gap end at the end of sequence

                                        std::cout<<"\t"<<pos+1<<"\n";

                                    }
                                    
                                    wasN = false;
                                    
                                    break;
                                    
                                }
                                    
                            }
                            
                            pos++;
                                    
                        }
                        
                        inSeq = "";
                        pos = 0;
                        gapLen = 0;
                        
                    }
                    
                }
                
                break;
                
            }
                
            case 'g': { // gaps
                
                unsigned int pos = 0, gapLen = 0;
                bool wasN = false;
                        
                std::string seqHeader, seqComment, inSeq; // header, comment and the new sequence being built recursively
                
                // generate adjacency list representation of a graph
                inSequences.buildGraph(inSequences.getGFAGaps());
                
                for (unsigned int i = 0; i != inSequences.getAdjListFW().size(); ++i) { // loop through all nodes
                    
                    InSegment inSegment; // a new inSequence object, the result of concatenating by gaps
                
                    seqHeader = inSequences.getInSegment(i).getSeqHeader();
                    seqComment = inSequences.getInSegment(i).getSeqComment();
                    
                    if (!inSequences.getVisited(i)) { // check if the node was already visited
                        
                        verbose(verbose_flag, "Graph DFS");
                        inSequences.dfsSeq(i, inSeq); // if not, visit all connected components recursively
                        
                        seqHeader.pop_back();
                        seqHeader.pop_back();
                        
                        for (char &base : inSeq) {

                            switch (base) {

                                case 'N': {
                                    
                                    gapLen++;
                                    
                                    if (!wasN) { // gap start
                                    
                                        std::cout<<seqHeader<<"\t"<<pos;
                                        
                                    }
                                    
                                    if (pos == inSeq.size()-1) { // gap end at the end of sequence
                                        
                                        std::cout<<"\t"<<pos+1<<"\n";
                                        
                                    }
                                    
                                    wasN = true;
                                    
                                    break;
                                    
                                }
                                default: {
                                    
                                    if (wasN) { // gap end
                                        
                                        if (pos == gapLen) { // gap at the start of sequence
                                        
                                            std::cout<<seqHeader<<"\t"<<0;
                                            
                                        }
                                    
                                        std::cout<<"\t"<<pos<<"\n";
                                        
                                    }
                                    
                                    wasN = false;
                                    
                                    break;
                                    
                                }
                                    
                            }
                            
                            pos++;
                                    
                        }
                        
                        inSeq = "";
                        pos = 0;
                        gapLen = 0;
                        
                    }
                    
                }
                
                break;
                
            }
                
            default:
            case 'a': { // both contigs and gaps in .agp format
                
                std::string seqHeader, seqComment, outAgp; // header, comment and the new sequence being built recursively
                unsigned int cStart = 1, cEnd = 1; // these are used to track coordinates along the scaffolds
  
                // generate adjacency list representation of a graph
                inSequences.buildGraph(inSequences.getGFAGaps());
                
                for (unsigned int i = 0; i != inSequences.getAdjListFW().size(); ++i) { // loop through all nodes
                    
                    InSegment inSegment; // a new inSequence object, the result of concatenating by gaps
                
                    seqHeader = inSequences.getInSegment(i).getSeqHeader();
                    seqHeader.pop_back();
                    seqHeader.pop_back();
                    
                    if (!inSequences.getVisited(i)) { // check if the node was already visited
                        
                        if (inSequences.getAdjListFW().at(i).size() == 0 && inSequences.getAdjListBW().at(i).size() == 0) { // handle disconnected components
                            
                            std::cout<<seqHeader<<"\t1\t"<<inSequences.getInSegment(i).getInSequence().size()<<"\t1\t"<<"W"<<"\t"<<seqHeader<<"\t1\t"<<inSequences.getInSegment(i).getInSequence().size()<<"\t+"<<"\n";
                            
                        }else{
                            
                            cStart = 1, cEnd = 1;
                            
                            verbose(verbose_flag, "Graph DFS");
                            
                            inSequences.dfsAgp(i, outAgp, cStart, cEnd); // if not, visit all connected components recursively
                        
                            std::cout<<outAgp;
                        
                        }
                    
                    }
                    
                    outAgp="";
                    
                }
                
                break;
            
            }
                
        }
        
        return true;
        
    }
    
    bool reportStats(InSequences &inSequences, unsigned long long int gSize, int bedOutType) { // method to output all summary statistic for the entire sequence set
        
        if (!tabular_flag) {
        
            std::cout<<output("+++Summary+++")<<"\n";
        
        }
        
        if (gSize > 0) {
        
            std::cout<<output("Expected genome size")<<gSize<<"\n";
        
        }
        
        std::cout<<output("# scaffolds")<<inSequences.getScaffN()<<"\n";
        std::cout<<output("Total scaffold length")<<inSequences.getTotScaffLen()<<"\n";
        printf("%s%.2f\n",output("Average scaffold length").c_str(), inSequences.computeAverageScaffLen());
        inSequences.evalNstars('s', gSize); // scaffold N* statistics
        std::cout<<output("Scaffold N50")<<inSequences.getScaffN50()<<"\n";
        inSequences.evalAuN('s', gSize); // scaffold auN
        printf("%s%.2f\n",output("Scaffold auN").c_str(), inSequences.getScaffauN());
        std::cout<<output("Scaffold L50")<<inSequences.getScaffL50()<<"\n";
        
        if (gSize > 0) {
            
            std::cout<<output("Scaffold NG50")<<inSequences.getScaffNG50()<<"\n";
            printf("%s%.2f\n",output("Scaffold auNG").c_str(), inSequences.getScaffauNG());
            std::cout<<output("Scaffold LG50")<<inSequences.getScaffLG50()<<"\n";
            
        }
        std::cout<<output("Largest scaffold")<<inSequences.getLargestScaffold()<<"\n";
        
        std::cout<<output("# contigs")<<inSequences.getSegmentN()<<"\n";
        std::cout<<output("Total contig length")<<inSequences.getTotSegmentLen()<<"\n";
        printf("%s%.2f\n",output("Average contig length").c_str(), inSequences.computeAverageSegmentLen());
        inSequences.evalNstars('c', gSize); // contig N* statistics
        std::cout<<output("Contig N50")<<inSequences.getContigN50()<<"\n";
        inSequences.evalAuN('c', gSize); // contig auN
        printf("%s%.2f\n",output("Contig auN").c_str(), inSequences.getContigauN());
        std::cout<<output("Contig L50")<<inSequences.getContigL50()<<"\n";
        
        if (gSize > 0) {
            
            std::cout<<output("Contig NG50")<<inSequences.getContigNG50()<<"\n";
            printf("%s%.2f\n",output("Contig auNG").c_str(), inSequences.getContigauNG());
            std::cout<<output("Contig LG50")<<inSequences.getContigLG50()<<"\n";
            
        }
        std::cout<<output("Largest contig")<<inSequences.getLargestContig()<<"\n";
        
        std::cout<<output("# gaps")<<inSequences.getGapN()<<"\n";
        std::cout<<output("Total gap length")<<inSequences.getTotGapLen()<<"\n";
        printf("%s%.2f\n",output("Average gap length").c_str(), inSequences.computeAverageGapLen());
        inSequences.evalNstars('g'); // gap N* statistics
        std::cout<<output("Gap N50")<<inSequences.getGapN50()<<"\n";
        inSequences.evalAuN('g'); // gap auN
        printf("%s%.2f\n",output("Gap auN").c_str(), inSequences.getGapauN());
        std::cout<<output("Gap L50")<<inSequences.getGapL50()<<"\n";
        std::cout<<output("Largest gap")<<inSequences.getLargestGap()<<"\n";
        
        printf("%s%lu, %lu, %lu, %lu\n",output("Base composition (ACGT)").c_str(), inSequences.getTotA(),
               inSequences.getTotC(),
               inSequences.getTotG(),
               inSequences.getTotT());
        printf("%s%.2f\n",output("GC content %").c_str(), inSequences.computeGCcontent());
        std::cout<<output("# soft-masked bases")<<inSequences.getTotLowerCount()<<"\n";
        
        counter = 0;
     
        return true;
        
    }
    
    bool nstarReport(InSequences &inSequences, unsigned long long int gSize) { // method to generate all N** reports
        
        int pos = 1;
        std::vector <unsigned int> scaffNstars = inSequences.getScaffNstars();
        for (unsigned int val : scaffNstars) {
            std::cout<<output("Scaffold N"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> scaffLstars = inSequences.getScaffLstars();
        for (unsigned int val : scaffLstars) {
            std::cout<<output("Scaffold L"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
        if (gSize > 0) {
            
            pos = 1;
            std::vector <unsigned int> scaffNGstars = inSequences.getScaffNGstars();
            for (unsigned int val : scaffNGstars) {
                std::cout<<output("Scaffold NG"+std::to_string(pos*10))<<val<<"\n";
                pos++;
            }
            
            pos = 1;
            std::vector <unsigned int> scaffLGstars = inSequences.getScaffLGstars();
            for (unsigned int val : scaffLGstars) {
                std::cout<<output("Scaffold LG"+std::to_string(pos*10))<<val<<"\n";
                pos++;
            }
            
        }
        
        pos = 1;
        std::vector <unsigned int> contigNstars = inSequences.getContigNstars();
        for (unsigned int val : contigNstars) {
            std::cout<<output("Contig N"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> contigLstars = inSequences.getContigLstars();
        for (unsigned int val : contigLstars) {
            std::cout<<output("Contig L"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
        if (gSize > 0) {
            
            pos = 1;
            std::vector <unsigned int> contigNGstars = inSequences.getContigNGstars();
            for (unsigned int val : contigNGstars) {
                std::cout<<output("Contig NG"+std::to_string(pos*10))<<val<<"\n";
                pos++;
            }
            
            pos = 1;
            std::vector <unsigned int> contigLGstars = inSequences.getContigLGstars();
            for (unsigned int val : contigLGstars) {
                std::cout<<output("Contig LG"+std::to_string(pos*10))<<val<<"\n";
                pos++;
            }
            
        }
        
        pos = 1;
        std::vector <unsigned int> gapNstars = inSequences.getGapNstars();
        for (unsigned int val : gapNstars) {
            std::cout<<output("Gap N"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> gapLstars = inSequences.getGapLstars();
        for (unsigned int val : gapLstars) {
            std::cout<<output("Gap L"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
        return true;
        
    }
    
    
};


#endif /* gfastats-output_h */
