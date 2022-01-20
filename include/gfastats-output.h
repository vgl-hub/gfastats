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
    bool seqReport(InSequences &inSequences, InSequence &inSequence, int &outSequence_flag) { // method to output the summary statistics for each sequence
        
        counter = 0;
        
        while (counter < inSequences.getScaffN()) {
            
            inSequence = inSequences.getInSequence(counter);
            
            std::cout<<output("Seq")<<counter+1<<std::endl;
            std::cout<<output("Header")<<inSequence.getSeqHeader()<<std::endl;
            std::cout<<output("Comment")<<inSequence.getSeqComment()<<std::endl;
            std::cout<<output("Total sequence length")<<inSequence.getSeqScaffLen()<<std::endl;
            std::cout<<output("Total contig length")<<inSequence.getContigSum()<<std::endl;
            std::cout<<output("# contig")<<inSequence.getContigN()<<std::endl;
            std::cout<<output("Total gap length")<<inSequence.getGapSum()<<std::endl;
            std::cout<<output("# gaps")<<inSequence.getGapN()<<std::endl;
            
            printf("%s%u, %u, %u, %u\n",output("Base composition (ACGT)").c_str(), inSequence.getA(),
                   inSequence.getC(),
                   inSequence.getG(),
                   inSequence.getT());
            printf("%s%.2f\n",output("GC content %").c_str(), inSequence.computeGCcontent());
            std::cout<<output("# soft-masked bases")<<inSequence.getLowerCount()<<std::endl;
            
            
            if (outSequence_flag) {
                
                std::cout<<output("Sequence")<<inSequence.getInSequence()<<std::endl;
                std::cout<<output("Quality")<<inSequence.getInSequenceQuality()<<std::endl;
                
            }
            
            std::cout<<std::endl;
            counter++;
            
        }
        
        counter = 0;
        
        return true;
        
    }
    
    bool outFile(InSequences &inSequences, InSequence &inSequence, int splitLength, std::string &outSeq) { // method to output new sequence opposed to sequence report
        
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
        
        // this stream outputs gzip compressed to stdout
        zstream::ogzstream zout(std::cout);
        
        // this stream outputs gzip to file
        std::ofstream ofs(outSeq);
        
        // this stream outputs gzip compressed to file
        zstream::ogzstream zfout(ofs);

        if (gzip && outFile) { // if the requested output is gzip compressed and should be outputted to a file
            
            stream = make_unique<std::ostream>(zfout.rdbuf()); // then we use the stream for gzip compressed file outputs
            
        }else if (!gzip && outFile){ // else if no compression is requested
            
            stream = make_unique<std::ostream>(ofs.rdbuf());  // we use the stream regular file outputs
            
        }else{ // else the output is not written to a file
            
            // we close and delete the file
            ofs.close();
            remove(outSeq.c_str());
            
            if (gzip) { // if the output to stdout needs to be compressed we use the appropriate stream
            
                stream = make_unique<std::ostream>(zout.rdbuf());
            
            }else{ // else we use a regular cout stream
                
                stream = make_unique<std::ostream>(std::cout.rdbuf());
                
            }
        }
        
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) { // this switch allows us to generate the output according to the input request and the unordered map. If the requested output format is not in the map we fall back to the undefined case (0)
                
            case 1: { // fasta[.gz]
                
                std::string output;
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSequence = inSequences.getInSequence(counter);
                    
                    *stream<<">"<<inSequence.getSeqHeader()<<" "<<inSequence.getSeqComment()<<std::endl;
                    
                    if (splitLength != 0) {
                        
                        textWrap(inSequence.getInSequence(), *stream, splitLength);
                        
                    }else{
                        
                        output = inSequence.getInSequence();
                        
                    }
                    
                    *stream<<output<<std::endl;
                    output = "";
                    
                    counter++;
                    
                }
                
                break;
                
            }
                
            case 2: { // fastq[.gz]
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSequence = inSequences.getInSequence(counter);
                    
                    *stream<<"@"<<inSequence.getSeqHeader()<<" "<<inSequence.getSeqComment()<<"\n"<<inSequence.getInSequence()<<"\n+\n"<<inSequence.getInSequenceQuality()<<"\n";
                    
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
                    
                    inSequence = inSequences.getInSequence(counter);
                    unsigned int seqScaffLen = inSequence.getSeqScaffLen();
                    
                    seqHeader = inSequence.getSeqHeader();
                    
                    seqBoundaries = inSequence.getSeqContigBoundaries();
                    
                    std::vector<unsigned int>::const_iterator begin = seqBoundaries.cbegin();
                    std::vector<unsigned int>::const_iterator end = seqBoundaries.cend();
                    auto last = std::prev(end);
                    
                    if (*begin>0) { // case apparently missing for GFA2 spec (starting gap)
                        
                        *stream <<"G\t" // line type
                                <<seqHeader<<"."<<item<<"\t" // id
                                <<seqHeader<<".begin"<<"\t" // sid1:ref (begin of sequence)
                                <<seqHeader<<"."<<item+1<<"\t" // sid2:ref
                                <<len<<"\t" // size
                                <<inSequence.getSeqComment()<<"\n"; // optional comment
                        
                        item++;
                        
                    }
                    
                    for (std::vector<unsigned int>::const_iterator it = seqBoundaries.cbegin(); it != end;) {
                        
                        len = *(it+1) - *it;

                        *stream <<"S\t" // line type
                                <<seqHeader<<"."<<item<<"\t" // id
                                <<*(it+1)-*(it)<<"\t" // size
                                <<inSequence.getInSequence().substr(*(it),*(it+1))<<"\t" // sequence
                                <<inSequence.getSeqComment()<<"\n"; // optional comment
                        
                        item++;
                        
                        if (ctgN != seqBoundaries.size()/2) {
                            
                            len = *(it+2) - *(it+1);
                            
                            *stream <<"G\t" // line type
                                    <<seqHeader<<"."<<item<<"\t" // id
                                    <<seqHeader<<"."<<item-1<<"\t" // sid1:ref
                                    <<seqHeader<<"."<<item+1<<"\t" // sid2:ref
                                    <<len<<"\t" // size
                                    <<inSequence.getSeqComment()<<"\n"; // optional comment
                            
                            item++;
                            
                        }
                        
                        if (ctgN == seqBoundaries.size()/2 && seqScaffLen > *last) {
                            
                            len = seqScaffLen - *(it+1);
                            
                            *stream <<"G\t" // line type
                                    <<seqHeader<<"."<<item<<"\t" // id
                                    <<seqHeader<<"."<<item-1<<"\t" // sid1:ref
                                    <<seqHeader<<".end"<<"\t" // sid2:ref (end of sequence)
                                    <<len<<"\t" // size
                                    <<inSequence.getSeqComment()<<"\n"; // optional comment
                            
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
    
    bool outSize(InSequences &inSequences, InSequence &inSequence, char &sizeOutType) { // method to output only the size of the sequences
        
        counter = 0;
        
        std::string seqHeader;
        std::vector<unsigned int> seqBoundaries;
        
        switch (sizeOutType) {
 
            default:
            case 's': { // scaffolds

                while (counter < inSequences.getScaffN()) {
                    
                    inSequence = inSequences.getInSequence(counter);
                        
                    std::cout<<inSequence.getSeqHeader()<<"\t"<<inSequence.getSeqScaffLen()<<std::endl;
                    
                    counter++;
                    
                }
                
                break;
            }
                
            case 'c': { // contigs
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSequence = inSequences.getInSequence(counter);
                    
                    seqHeader = inSequence.getSeqHeader();
                    
                    seqBoundaries = inSequence.getSeqContigBoundaries();
                    
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
                    
                    inSequence = inSequences.getInSequence(counter);
                    
                    seqHeader = inSequence.getSeqHeader();
                    
                    seqBoundaries = inSequence.getSeqGapBoundaries();
                    
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
    
    bool outCoord(InSequences &inSequences, InSequence &inSequence, char bedOutType) { // method to output the coordinates of each feature
        
        counter = 0;
        
        std::string seqHeader;
        std::vector<unsigned int> seqBoundaries;
        
        switch (bedOutType) {
                
            case 'c': { // contigs
                
                while (counter < inSequences.getScaffN()) {
                    
                    inSequence = inSequences.getInSequence(counter);
                    
                    seqHeader = inSequence.getSeqHeader();
                    
                    seqBoundaries = inSequence.getSeqContigBoundaries();
                    
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
                    
                    inSequence = inSequences.getInSequence(counter);
                    
                    seqHeader = inSequence.getSeqHeader();
                    
                    seqBoundaries = inSequence.getSeqGapBoundaries();
                    
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
                    
                    inSequence = inSequences.getInSequence(counter);
                    unsigned int seqScaffLen = inSequence.getSeqScaffLen();
                    
                    seqHeader = inSequence.getSeqHeader();
                    
                    seqBoundaries = inSequence.getSeqContigBoundaries();
                    
                    std::vector<unsigned int>::const_iterator begin = seqBoundaries.cbegin();
                    std::vector<unsigned int>::const_iterator end = seqBoundaries.cend();
                    auto last = std::prev(end);
                    
                    if (*begin>0) {
                        
                        std::cout<<seqHeader<<"\t"<<1<<"\t"<<*begin<<"\t"<<1<<"\t"<<"N"<<"\t"<<*begin<<"\tscaffold\tyes\t"<<std::endl;
                        
                        item++;
                        
                    }
                    
                    for (std::vector<unsigned int>::const_iterator it = seqBoundaries.cbegin(); it != end;) {
                        
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
        std::cout<<output("Total contig length")<<inSequences.getTotContigLen()<<std::endl;
        printf("%s%.2f\n",output("Average contig length").c_str(), inSequences.computeAverageContigLen());
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
