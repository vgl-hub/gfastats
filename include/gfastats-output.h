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
        
//        counter = 0;
//
//        while (counter < inSequences.getScaffN()) {
//
//            inSegment = inSequences.getInSequence(counter);
//
//            std::cout<<output("Seq")<<counter+1<<std::endl;
//            std::cout<<output("Header")<<inSegment.getSeqHeader()<<std::endl;
//            std::cout<<output("Comment")<<inSegment.getSeqComment()<<std::endl;
//            std::cout<<output("Total sequence length")<<inSegment.getSegmentLength()<<std::endl;
//            std::cout<<output("Total contig length")<<inSegment.getContigSum()<<std::endl;
//            std::cout<<output("# contig")<<inSegment.getContigN()<<std::endl;
//            std::cout<<output("Total gap length")<<inSegment.getGapSum()<<std::endl;
//            std::cout<<output("# gaps")<<inSegment.getGapN()<<std::endl;
//
//            printf("%s%u, %u, %u, %u\n",output("Base composition (ACGT)").c_str(), inSequence.getA(),
//                   inSequence.getC(),
//                   inSequence.getG(),
//                   inSequence.getT());
//            printf("%s%.2f\n",output("GC content %").c_str(), inSequence.computeGCcontent());
//            std::cout<<output("# soft-masked bases")<<inSequence.getLowerCount()<<std::endl;
//
//
//            if (outSequence_flag) {
//
//                std::cout<<output("Sequence")<<inSequence.getInSequence()<<std::endl;
//                std::cout<<output("Quality")<<inSequence.getInSequenceQuality()<<std::endl;
//
//            }
//
//            std::cout<<std::endl;
//            counter++;
//
//        }
//
//        counter = 0;
        
        return true;
        
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
                
                std::string output;
                
                std::string seqHeader, seqComment, inSeq;
                
                // gfa postprocessing
                // print adjacency list representation of a graph
                inSequences.buildGraph(inSequences.getGFAGaps());
                
                InSequences inSequencesNew; // the new sequence set resulting from the graph
                
                
                for (unsigned int i = 0; i != inSequences.getAdjListFW().size(); ++i) { // loop through all nodes
                    
                    InSegment inSegment; // a new inSequence object, the result of concatenating by gaps
                    std::string inSeqNew; // the new sequence being built recursively
                
                    seqHeader = inSequences.getInSegment(i).getSeqHeader();
                    seqComment = inSequences.getInSegment(i).getSeqComment();
                    
                    if (!inSequences.getVisited(i)) { // check if the node was already visited
                        
                        inSequences.DFS(i, inSegment, inSeqNew); // if not, visit all connected components recursively
                        
                        inSeq = inSegment.getInSegment();
                    
                        inSequencesNew.addSegment(&seqHeader, &seqComment, &inSeq); // push the new sequence
                        
                        inSequencesNew.increaseTotScaffN();
                        
                    }
                    
                }
                
                inSequences = inSequencesNew;
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                    
                    *stream<<">"<<inSegment.getSeqHeader()<<" "<<inSegment.getSeqComment()<<std::endl;
                    
                    if (splitLength != 0) {
                        
                        textWrap(inSegment.getInSegment(), *stream, splitLength);
                        
                    }else{
                        
                        output = inSegment.getInSegment();
                        
                    }
                    
                    *stream<<output<<std::endl;
                    output = "";
                    
                    counter++;
                    
                }
                
                break;
                
            }
                
            case 2: { // fastq[.gz]
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                    
                    *stream<<"@"<<inSegment.getSeqHeader()<<" "<<inSegment.getSeqComment()<<"\n"<<inSegment.getInSegment()<<"\n+\n"<<inSegment.getInSequenceQuality()<<"\n";
                    
                    counter++;
                    
                }
                
                break;
                
            }
                
            case 3: { // gfa[.gz]
                
                std::string seqHeader;
                std::vector<unsigned int> seqBoundaries;
                unsigned int ctgN = 1, item = 1, len = 0;
                
                *stream<<"H\tVN:Z:2.0\n";
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                    unsigned int seqScaffLen = inSegment.getSegmentLength();
                    
                    seqHeader = inSegment.getSeqHeader();
                    
                    std::vector<unsigned int>::const_iterator begin = seqBoundaries.cbegin();
                    std::vector<unsigned int>::const_iterator end = seqBoundaries.cend();
                    auto last = std::prev(end);
                    
                    if (*begin>0) { // case apparently missing for GFA2 spec (start gap)
                        
                        len = *begin;
                        
                        *stream <<"G\t" // line type
                                <<seqHeader<<"."<<item<<"\t" // id
                                <<seqHeader<<"."<<item+1<<"+\t" // sid1:ref (begin of sequence)
                                <<seqHeader<<"."<<item+1<<"+\t" // sid2:ref
                                <<len<<"\t" // size
                                <<inSegment.getSeqComment()<<"\n"; // optional comment
                        
                        item++;
                        
                    }
                    
                    for (std::vector<unsigned int>::const_iterator it = seqBoundaries.cbegin(); it != end;) {
                        
                        len = *(it+1) - *it;

                        *stream <<"S\t" // line type
                                <<seqHeader<<"."<<item<<"\t" // id
                                <<*(it+1)-*(it)<<"\t" // size
                                <<inSegment.getInSegment().substr(*(it),len)<<"\t" // sequence
                                <<inSegment.getSeqHeader(); // header
                        
                        if (inSegment.getSeqComment() != "") {
                        
                        *stream <<" "<<inSegment.getSeqComment(); // optional comment
                    
                        }
                
                        *stream <<"\n";
                        
                        item++;
                        
                        if (ctgN != seqBoundaries.size()/2) { // inner gap
                            
                            len = *(it+2) - *(it+1);
                            
                            *stream <<"G\t" // line type
                                    <<seqHeader<<"."<<item<<"\t" // id
                                    <<seqHeader<<"."<<item-1<<"+\t" // sid1:ref
                                    <<seqHeader<<"."<<item+1<<"+\t" // sid2:ref
                                    <<len<<"\t" // size
                                    <<inSegment.getSeqComment()<<"\n"; // optional comment
                            
                            item++;
                            
                        }
                        
                        if (ctgN == seqBoundaries.size()/2 && seqScaffLen > *last) { // end gap
                            
                            len = seqScaffLen - *(it+1);
                            
                            *stream <<"G\t" // line type
                                    <<seqHeader<<"."<<item<<"\t" // id
                                    <<seqHeader<<"."<<item-1<<"+\t" // sid1:ref
                                    <<seqHeader<<"."<<item-1<<"-\t" // sid2:ref (end of sequence)
                                    <<len<<"\t" // size
                                    <<inSegment.getSeqComment()<<"\n"; // optional comment
                            
                            item++;
                            
                        }
                        
                        ctgN++;
                        it = it + 2;
                        
                    }
                    
                    ctgN = 1;
                    item = 1;
                    counter++;
                    
                }
                
                counter = 0;
                std::vector<InGap> inGaps = inSequences.getGFAGaps();
                
                while (counter < inGaps.size()) {
                
                    *stream <<"G\t" // line type
                            <<inGaps[counter].getgIds()<<"\t" // id
                            <<inGaps[counter].getsIds1()<<"\t" // sid1:ref
                            <<inGaps[counter].getsIds2()<<"\t" // sid2:ref (end of sequence)
                            <<inGaps[counter].getsDists()<<"\n"; // size
                    
                    counter++;
                    
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
        
        counter = 0;
        
        std::string seqHeader;
        std::vector<unsigned int> seqBoundaries;
        
        switch (sizeOutType) {
 
            default:
            case 's': { // scaffolds

                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                        
                    std::cout<<inSegment.getSeqHeader()<<"\t"<<inSegment.getSegmentLength()<<std::endl;
                    
                    counter++;
                    
                }
                
                break;
            }
                
            case 'c': { // contigs
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                    
                    seqHeader = inSegment.getSeqHeader();
                    
                    std::vector<unsigned int>::const_iterator end = seqBoundaries.cend();
                    
                    for (std::vector<unsigned int>::const_iterator it = seqBoundaries.cbegin(); it != end;) {
                        
                        std::cout<<seqHeader<<"\t"<<*(it+1)-*it<<std::endl;
                        
                        it = it + 2;
                        
                    }
                    
                    counter++;
                    
                }
                
                break;
                
            }
                
            case 'g': { // gaps
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                    
                    seqHeader = inSegment.getSeqHeader();
                    
                    std::vector<unsigned int>::const_iterator end = seqBoundaries.cend();
                    
                    for (std::vector<unsigned int>::const_iterator it = seqBoundaries.cbegin(); it != end;) {
                        
                        std::cout<<seqHeader<<"\t"<<*(it+1)-*it<<std::endl;
                        
                        it = it + 2;
                        
                    }
                    
                    counter++;
                    
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
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                    
                    seqHeader = inSegment.getSeqHeader();
                    
                    std::vector<unsigned int>::const_iterator end = seqBoundaries.cend();
                    
                    for (std::vector<unsigned int>::const_iterator it = seqBoundaries.cbegin(); it != end;) {
                        
                        std::cout<<seqHeader<<"\t"<<*it<<"\t"<<*(it+1)<<std::endl;
                        
                        it = it + 2;
                        
                    }
                    
                    counter++;
                    
                }
                
                break;
                
            }
                
            case 'g': { // gaps
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                    
                    seqHeader = inSegment.getSeqHeader();
                    
                    std::vector<unsigned int>::const_iterator end = seqBoundaries.cend();
                    
                    for (std::vector<unsigned int>::const_iterator it = seqBoundaries.cbegin(); it != end;) {
                        
                        std::cout<<seqHeader<<"\t"<<*it<<"\t"<<*(it+1)<<std::endl;
                        
                        it = it + 2;
                        
                    }
                    
                    counter++;
                    
                }
                
                break;
                
            }
                
            default:
            case 'a': { // both contigs and gaps in .agp format
                
                unsigned int ctgN = 1, item = 1, len = 0;
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSegment = inSequences.getInSegment(counter);
                    unsigned int seqScaffLen = inSegment.getSegmentLength();
                    
                    seqHeader = inSegment.getSeqHeader();
                    
                    std::vector<unsigned int>::const_iterator begin = seqBoundaries.cbegin();
                    std::vector<unsigned int>::const_iterator end = seqBoundaries.cend();
                    auto last = std::prev(end);
                    
                    if (*begin>0) {
                        
                        std::cout<<seqHeader<<"\t"<<1<<"\t"<<*begin<<"\t"<<1<<"\t"<<"N"<<"\t"<<*begin<<"\tscaffold\tyes\t"<<std::endl;
                        
                        item++;
                        
                    }
                    
                    for (std::vector<unsigned int>::const_iterator it = begin; it != end;) {
                        
                        len = *(it+1) - *it;
                        
                        std::cout<<seqHeader<<"\t"<<*it+1<<"\t"<<*(it+1)<<"\t"<<item<<"\t"<<"W"<<"\t"<<seqHeader+"."<<ctgN<<"\t1\t"<<len<<"\t+"<<std::endl;
                        
                        item++;
                        
                        if (ctgN != seqBoundaries.size()/2) {
                            
                            len = *(it+2) - *(it+1);
                            
                            std::cout<<seqHeader<<"\t"<<*(it+1)+1<<"\t"<<*(it+2)<<"\t"<<item<<"\t"<<"N"<<"\t"<<len<<"\tscaffold\tyes\t"<<std::endl;
                            
                            item++;
                            
                        }
                        
                        if (ctgN == seqBoundaries.size()/2 && seqScaffLen > *last) {
                            
                            len = seqScaffLen - *(it+1);
                            
                            std::cout<<seqHeader<<"\t"<<*(it+1)+1<<"\t"<<seqScaffLen<<"\t"<<item<<"\t"<<"N"<<"\t"<<len<<"\tscaffold\tyes\t"<<std::endl;
                            
                            item++;
                            
                        }
                        
                        ctgN++;
                        it = it + 2;
                        
                    }
                    
                    ctgN = 1;
                    item = 1;
                    counter++;
                    
                }
                
                break;
            }
                
        }
        
        return true;
        
    }
    
    bool reportStats(InSequences &inSequences, unsigned long long int gSize, int bedOutType) { // method to output all summary statistic for the entire sequence set
        
        if (!tabular_flag) {
        
            std::cout<<output("+++Summary+++")<<std::endl;
        
        }
        
        if (gSize > 0) {
        
            std::cout<<output("Expected genome size")<<gSize<<std::endl;
        
        }
        
        std::cout<<output("# scaffolds")<<inSequences.getScaffN()<<std::endl;
        std::cout<<output("Total scaffold length")<<inSequences.getTotScaffLen()<<std::endl;
        printf("%s%.2f\n",output("Average scaffold length").c_str(), inSequences.computeAverageScaffLen());
        inSequences.evalNstars('s', gSize); // scaffold N* statistics
        std::cout<<output("Scaffold N50")<<inSequences.getScaffN50()<<std::endl;
        inSequences.evalAuN('s', gSize); // scaffold auN
        printf("%s%.2f\n",output("Scaffold auN").c_str(), inSequences.getScaffauN());
        std::cout<<output("Scaffold L50")<<inSequences.getScaffL50()<<std::endl;
        
        if (gSize > 0) {
            
            std::cout<<output("Scaffold NG50")<<inSequences.getScaffNG50()<<std::endl;
            printf("%s%.2f\n",output("Scaffold auNG").c_str(), inSequences.getScaffauNG());
            std::cout<<output("Scaffold LG50")<<inSequences.getScaffLG50()<<std::endl;
            
        }
        std::cout<<output("Largest scaffold")<<inSequences.getLargestScaffold()<<std::endl;
        
        std::cout<<output("# contigs")<<inSequences.getContigN()<<std::endl;
        std::cout<<output("Total contig length")<<inSequences.getTotSegmentLen()<<std::endl;
        printf("%s%.2f\n",output("Average contig length").c_str(), inSequences.computeAverageSegmentLen());
        inSequences.evalNstars('c', gSize); // contig N* statistics
        std::cout<<output("Contig N50")<<inSequences.getContigN50()<<std::endl;
        inSequences.evalAuN('c', gSize); // contig auN
        printf("%s%.2f\n",output("Contig auN").c_str(), inSequences.getContigauN());
        std::cout<<output("Contig L50")<<inSequences.getContigL50()<<std::endl;
        
        if (gSize > 0) {
            
            std::cout<<output("Contig NG50")<<inSequences.getContigNG50()<<std::endl;
            printf("%s%.2f\n",output("Contig auNG").c_str(), inSequences.getContigauNG());
            std::cout<<output("Contig LG50")<<inSequences.getContigLG50()<<std::endl;
            
        }
        std::cout<<output("Largest contig")<<inSequences.getLargestContig()<<std::endl;
        
        std::cout<<output("# gaps")<<inSequences.getTotGapN()<<std::endl;
        std::cout<<output("Total gap length")<<inSequences.getTotGapLen()<<std::endl;
        printf("%s%.2f\n",output("Average gap length").c_str(), inSequences.computeAverageGapLen());
        inSequences.evalNstars('g'); // gap N* statistics
        std::cout<<output("Gap N50")<<inSequences.getGapN50()<<std::endl;
        inSequences.evalAuN('g'); // gap auN
        printf("%s%.2f\n",output("Gap auN").c_str(), inSequences.getGapauN());
        std::cout<<output("Gap L50")<<inSequences.getGapL50()<<std::endl;
        std::cout<<output("Largest gap")<<inSequences.getLargestGap()<<std::endl;
        
        printf("%s%lu, %lu, %lu, %lu\n",output("Base composition (ACGT)").c_str(), inSequences.getTotA(),
               inSequences.getTotC(),
               inSequences.getTotG(),
               inSequences.getTotT());
        printf("%s%.2f\n",output("GC content %").c_str(), inSequences.computeGCcontent());
        std::cout<<output("# soft-masked bases")<<inSequences.getTotLowerCount()<<std::endl;
        
        counter = 0;
     
        return true;
        
    }
    
    bool nstarReport(InSequences &inSequences, unsigned long long int gSize) { // method to generate all N** reports
        
        int pos = 1;
        std::vector <unsigned int> scaffNstars = inSequences.getScaffNstars();
        for (unsigned int val : scaffNstars) {
            std::cout<<output("Scaffold N"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> scaffLstars = inSequences.getScaffLstars();
        for (unsigned int val : scaffLstars) {
            std::cout<<output("Scaffold L"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        if (gSize > 0) {
            
            pos = 1;
            std::vector <unsigned int> scaffNGstars = inSequences.getScaffNGstars();
            for (unsigned int val : scaffNGstars) {
                std::cout<<output("Scaffold NG"+std::to_string(pos*10))<<val<<std::endl;
                pos++;
            }
            
            pos = 1;
            std::vector <unsigned int> scaffLGstars = inSequences.getScaffLGstars();
            for (unsigned int val : scaffLGstars) {
                std::cout<<output("Scaffold LG"+std::to_string(pos*10))<<val<<std::endl;
                pos++;
            }
            
        }
        
        pos = 1;
        std::vector <unsigned int> contigNstars = inSequences.getContigNstars();
        for (unsigned int val : contigNstars) {
            std::cout<<output("Contig N"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> contigLstars = inSequences.getContigLstars();
        for (unsigned int val : contigLstars) {
            std::cout<<output("Contig L"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        if (gSize > 0) {
            
            pos = 1;
            std::vector <unsigned int> contigNGstars = inSequences.getContigNGstars();
            for (unsigned int val : contigNGstars) {
                std::cout<<output("Contig NG"+std::to_string(pos*10))<<val<<std::endl;
                pos++;
            }
            
            pos = 1;
            std::vector <unsigned int> contigLGstars = inSequences.getContigLGstars();
            for (unsigned int val : contigLGstars) {
                std::cout<<output("Contig LG"+std::to_string(pos*10))<<val<<std::endl;
                pos++;
            }
            
        }
        
        pos = 1;
        std::vector <unsigned int> gapNstars = inSequences.getGapNstars();
        for (unsigned int val : gapNstars) {
            std::cout<<output("Gap N"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> gapLstars = inSequences.getGapLstars();
        for (unsigned int val : gapLstars) {
            std::cout<<output("Gap L"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        return true;
        
    }
    
    
};


#endif /* gfastats-output_h */
