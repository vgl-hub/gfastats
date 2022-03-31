//
//  gfastats-output.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFASTATS_OUTPUT_H
#define GFASTATS_OUTPUT_H

//classes
class Report {
private:
    unsigned int counter = 0;
    
public:
    bool seqReport(InSequences &inSequences, InSegment &inSegment, int &outSequence_flag) { // method to output the summary statistics for each sequence
        std::cout << std::fixed; // disables scientific notation
        std::cout << std::setprecision(2); // 2 decimal poinst
        counter = 0;
        std::vector<InSegment>* inSegments = inSequences.getInSegments();
        
        std::cout<<output("Seq\tHeader\tComment\tTotal segment length\tA\tC\tG\tT\tGC content %\t# soft-masked bases");
        
        if (outSequence_flag) {

            std::cout<<output("Sequence\tQuality");

        }
        
        for (InSegment inSegment : *inSegments) {
            
            std::cout   <<"\n"<<counter+1<<"\t"
                        <<inSegment.getSeqHeader()<<"\t"
                        <<inSegment.getSeqComment()<<"\t"
                        <<inSegment.getSegmentLen()<<"\t"
                        <<inSegment.getA()<<"\t"
                        <<inSegment.getC()<<"\t"
                        <<inSegment.getG()<<"\t"
                        <<inSegment.getT()<<"\t"
                        <<inSegment.computeGCcontent()<<"\t"
                        <<inSegment.getLowerCount();

            if (outSequence_flag) {

                std::cout<<inSegment.getInSequence()<<"\t"<<inSegment.getInSequenceQuality()<<"\n";

            }
            
            counter++;

        }
        
        std::cout<<"\n";

        counter = 0;
        
        return true;
        
    }
    
    bool outFile(InSequences &inSequences, InSegment &inSegment, int splitLength, std::string &outSeq) { // method to output new sequence opposed to sequence report
        std::cout << std::fixed; // disables scientific notation
        std::cout << std::setprecision(2); // 2 decimal poinst

        // unordered map to handle out correspondence in following switch statement
        const static phmap::flat_hash_map<std::string,int> string_to_case{
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
            
            stats_flag = true; // since we write to file, let's output the stats
            
        }else{
            
            ext = outSeq;
            
        }
        
        // here we create a smart pointer to handle any kind of output stream
        std::unique_ptr<std::ostream> stream;
        
        // this stream outputs to file
        std::ofstream ofs(outSeq);

        // this stream outputs gzip compressed to file
        zstream::ogzstream zfout(ofs);
        
        // this stream outputs gzip compressed to stdout
        zstream::ogzstream zout(std::cout);

        if (gzip && outFile) { // if the requested output is gzip compressed and should be outputted to a file
            
            stream = make_unique<std::ostream>(zfout.rdbuf()); // then we use the stream for gzip compressed file outputs
            
            zfout.addHeader();
                        
        }else if (!gzip && outFile){ // else if no compression is requested
            
            stream = make_unique<std::ostream>(ofs.rdbuf());  // we use the stream regular file outputs
            
        }else{ // else the output is not written to a file
            
            // we close and delete the file
            ofs.close();
            remove(outSeq.c_str());
            
            if (gzip) { // if the output to stdout needs to be compressed we use the appropriate stream
                
                stream = make_unique<std::ostream>(zout.rdbuf());
                
                zout.addHeader();
            
            }else{ // else we use a regular cout stream
                
                std::cout.flush();
                
                stream = make_unique<std::ostream>(std::cout.rdbuf());
                
            }
        }
        
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) { // this switch allows us to generate the output according to the input request and the unordered map. If the requested output format is not in the map we fall back to the undefined case (0)
                
            case 1: { // fasta[.gz]
                
                std::string pHeader;
                std::string inSeq; // the new sequence being built
                std::vector<InPath> inPaths = inSequences.getInPaths();
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                std::vector<InGap>* inGaps = inSequences.getInGaps();
                std::vector<PathTuple> pathComponents;
                
                unsigned int uId = 0, sIdx = 0, gIdx = 0;
                    
                for (InPath inPath : inSequences.getInPaths()) {
                    
                    if (inPath.getHeader() == "") {
                        
                        pHeader = inPath.getpUId();
                        
                    }else{
                        
                        pHeader = inPath.getHeader();
                        
                    }
                    
                    *stream <<">"
                            <<pHeader;
                    
                    if (inPath.getComment() != "") {
                    
                    *stream <<" "
                            <<inPath.getComment();
                        
                    }
                    
                    *stream <<"\n";
                    
                    pathComponents = inPath.getComponents();
                    
                    for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                        
                        uId = std::get<1>(*component);
                        
                        if (std::get<0>(*component) == 'S') {
                        
                            auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                            
                            if (std::get<2>(*component) == '+') {
                            
                                inSeq += (*inSegments)[sIdx].getInSequence(std::get<3>(*component), std::get<4>(*component));
                                
                            }else{
                                
                                inSeq += revCom((*inSegments)[sIdx].getInSequence(std::get<3>(*component), std::get<4>(*component)));
                                
                            }
                            
                        }else{
                            
                            auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the segment index
                            
                            inSeq += std::string((*inGaps)[gIdx].getDist(), 'N');
                            
                        }
                        
                    }
                    
                    if (splitLength != 0) {
                        
                        textWrap(inSeq, *stream, splitLength); // wrapping text at user-specified line length
                        
                    }else{
                        
                        *stream<<inSeq<<"\n";
                        
                    }
                    
                    inSeq = "";
                        
                    (*stream).flush();
                    
                }
                
                break;
                
            }
                
            case 2: { // fastq[.gz]
                
                std::string pHeader;
                std::string inSeq, inSeqQual; // the new sequence being built and its quality
                std::vector<InPath> inPaths = inSequences.getInPaths();
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                std::vector<InGap>* inGaps = inSequences.getInGaps();
                std::vector<PathTuple> pathComponents;
                
                unsigned int uId = 0, sIdx = 0, gIdx = 0;
                    
                for (InPath inPath : inSequences.getInPaths()) {
                    
                    if (inPath.getHeader() == "") {
                        
                        pHeader = inPath.getpUId();
                        
                    }else{
                        
                        pHeader = inPath.getHeader();
                        
                    }
                    
                    *stream <<"@"
                            <<pHeader;
                    
                    if (inPath.getComment() != "") {
                    
                    *stream <<" "
                            <<inPath.getComment();
                        
                    }
                    
                    *stream <<"\n";
                    
                    pathComponents = inPath.getComponents();
                    
                    for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                        
                        uId = std::get<1>(*component);
                        
                        if (std::get<0>(*component) == 'S') {
                        
                            auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                            
                            if (std::get<2>(*component) == '+') {
                            
                                inSeq += (*inSegments)[sIdx].getInSequence(std::get<3>(*component), std::get<4>(*component));
                                
                            }else{
                                
                                inSeq += revCom((*inSegments)[sIdx].getInSequence(std::get<3>(*component), std::get<4>(*component)));
                                
                            }
                            
                            if ((*inSegments)[sIdx].getInSequenceQuality() != "") {
                            
                                inSeqQual += (*inSegments)[sIdx].getInSequence(std::get<3>(*component), std::get<4>(*component));
                                
                            }else{
                                
                                inSeqQual += std::string((*inSegments)[sIdx].getInSequence().size(), '!');
                                
                            }
                            
                        }else{
                            
                            auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the segment index
                            
                            inSeq += std::string((*inGaps)[gIdx].getDist(), 'N');
                            inSeqQual += std::string((*inGaps)[gIdx].getDist(), '!');
                            
                        }
                        
                        
                    }
                    
                    *stream<<inSeq<<"\n+\n"<<inSeqQual<<"\n";
                    
                    inSeq = "";
                    inSeqQual = "";
                    
                }
                
                break;
                
            }
                
            case 3: { // gfa[.gz]
                
                std::string seqHeader, gHeader, pHeader;
                
                phmap::flat_hash_map<unsigned int, std::string> idsToHeaders = inSequences.getHash2();
                
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                
                *stream<<"H\tVN:Z:2.0\n";
                
                for (InSegment inSegment : *inSegments) {
                    
                    seqHeader = inSegment.getSeqHeader();
                    
                    *stream <<"S\t" // line type
                            <<seqHeader<<"\t" // header
                            <<inSegment.getSegmentLen()<<"\t" // seq length
                            <<inSegment.getInSequence(); // sequence
                    
                    if (inSegment.getInSequenceQuality() != "") {
                        
                        *stream <<"\tQ:"<<inSegment.getInSequenceQuality(); // optional quality
                        
                    }
                    
                    *stream<<"\n";
                    
                }
                
                for (InGap inGap : inSequences.getGaps()) {
                    
                    if (inGap.getgHeader() == "") {
                        
                        gHeader = inGap.getuId();
                        
                    }else{
                        
                        gHeader = inGap.getgHeader();
                        
                    }
                    
                    *stream <<"G\t" // line type
                            <<gHeader<<"\t" // id
                            <<idsToHeaders[inGap.getsId1()]<<inGap.getsId1Or()<<"\t" // sUid1:sid1:ref
                            <<idsToHeaders[inGap.getsId2()]<<inGap.getsId2Or()<<"\t" // sUid2:sid2:ref
                            <<inGap.getDist()<<"\n"; // size
                    
                }
                
                std::vector<PathTuple> pathComponents;
                
                for (InPath inPath : inSequences.getInPaths()) {
                    
                    if (inPath.getHeader() == "") {
                        
                        pHeader = inPath.getpUId();
                        
                    }else{
                        
                        pHeader = inPath.getHeader();
                        
                    }
                    
                    *stream <<"O\t" // line type
                            <<pHeader<<"\t"; // id
                    
                    
                    pathComponents = inPath.getComponents();
                    
                    for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                            
                        *stream << idsToHeaders[std::get<1>(*component)];
                        
                        if(std::get<3>(*component) != 0) {
                        
                            *stream << "(" << std::to_string(std::get<3>(*component)) << ":" << std::to_string(std::get<4>(*component)) << ")";
                            
                        }
                        
                        if(std::get<2>(*component) != '0') {
                        
                            *stream << std::get<2>(*component);
                            
                        }
                        
                        if (component != std::prev(pathComponents.end())) {
                            
                            *stream <<" "; // space
                            
                        }
                        
                    }
                    
                    if (inPath.getComment() != "") {
                    
                    *stream <<"\t"
                            <<inPath.getComment();
                        
                    }
                    
                    *stream <<"\n";
                    
                }
                
                break;
                
            }
                
            case 0: { // undefined case
                
                std::cout<<"Unrecognized output format: "<<outSeq;
                
                break;
                
            }
                
        }
        
        if(gzip && outFile) { // if we wrote to file as gzip, we add the footer and close
            
            zfout.close();
            
        }else if(gzip && !outFile) { // if we streamed as gzip, we add the footer and close
            
            zout.close();
            
        }
        
        if(outFile) { // if we wrote to file, we close it
            
            ofs.close();
            
        }
        
        return true;
        
    }
    
    bool outSize(InSequences &inSequences, InSegment &inSegment, char &sizeOutType) { // method to output only the size of the sequences
        
        std::cout << std::fixed; // disables scientific notation
        std::cout << std::setprecision(2); // 2 decimal poinst
        switch (sizeOutType) {
 
            default:
            case 's': { // scaffolds
                
                std::string pHeader;
                std::vector<InPath> inPaths = inSequences.getInPaths();
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                std::vector<InGap>* inGaps = inSequences.getInGaps();
                std::vector<PathTuple> pathComponents;
                
                unsigned int uId = 0, sIdx = 0, gIdx = 0, size = 0;
                    
                for (InPath inPath : inSequences.getInPaths()) {
                    
                    if (inPath.getHeader() == "") {
                        
                        pHeader = inPath.getpUId();
                        
                    }else{
                        
                        pHeader = inPath.getHeader();
                        
                    }
                    
                    std::cout<<pHeader<<"\t";
                    
                    pathComponents = inPath.getComponents();
                    
                    for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                        
                        if (std::get<4>(*component) > 0) { // if we are subsetting we do not need to know the length of the segment
                            
                            size += std::get<4>(*component) - std::get<3>(*component);
                            
                        }else{
                        
                            uId = std::get<1>(*component);
                            
                            if (std::get<0>(*component) == 'S') {
                            
                                auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                                
                                if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                                
                                size += (*inSegments)[sIdx].getInSequence().size();
                                
                            }else{
                                
                                auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                                
                                if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the segment index
                                
                                size += (*inGaps)[gIdx].getDist();
                                
                            }
                            
                        }
                        
                    }
                    
                    std::cout<<size<<"\n";
                    size = 0;
                    
                }
                
                break;
            }
                
            case 'c': { // contigs
                
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                
                for (InSegment inSegment : *inSegments) {
                    
                    std::cout<<inSegment.getSeqHeader()<<"\t"<<inSegment.getInSequence().size()<<"\n";

                }
                
                break;
                
            }
                
            case 'g': { // gaps
                
                std::vector<InGap>* inGaps = inSequences.getInGaps();
                
                for (InGap inGap : *inGaps) {

                    std::cout<<inGap.getgHeader()<<"\t"<<inGap.getDist()<<"\n";

                }
                
                break;
                
            }
                
        }
        
        return true;
        
    }
    
    bool outCoord(InSequences &inSequences, InSegment &inSegment, char bedOutType) { // method to output the coordinates of each feature
        std::cout << std::fixed; // disables scientific notation
        std::cout << std::setprecision(2); // 2 decimal poinst

        counter = 0;
        
        std::string seqHeader;
        std::vector<unsigned int> seqBoundaries;
        
        switch (bedOutType) {
                
            case 's': { // scaffolds
                
                std::string pHeader;
                std::vector<InPath> inPaths = inSequences.getInPaths();
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                std::vector<InGap>* inGaps = inSequences.getInGaps();
                std::vector<PathTuple> pathComponents;
                
                unsigned int uId = 0, sIdx = 0, gIdx = 0, pos = 0;
                    
                for (InPath inPath : inSequences.getInPaths()) {
                    
                    if (inPath.getHeader() == "") {
                        
                        pHeader = inPath.getpUId();
                        
                    }else{
                        
                        pHeader = inPath.getHeader();
                        
                    }
                    
                    std::cout<<pHeader<<"\t"<<pos;
                    
                    pathComponents = inPath.getComponents();
                    
                    for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                        
                        uId = std::get<1>(*component);
                        
                        if (std::get<0>(*component) == 'S') {
                        
                            auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                            
                            pos += (*inSegments)[sIdx].getInSequence().size();
                            
                        }else{
                            
                            auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the segment index
                            
                            pos += (*inGaps)[gIdx].getDist();
                            
                        }
                        
                        
                    }
                    
                    std::cout<<"\t"<<pos<<"\n";
                    pos = 0;
                    
                }
                
                break;
                
            }
                
            case 'c': { // contigs
                
                std::string pHeader;
                std::vector<InPath> inPaths = inSequences.getInPaths();
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                std::vector<InGap>* inGaps = inSequences.getInGaps();
                std::vector<PathTuple> pathComponents;
                
                unsigned int uId = 0, sIdx = 0, gIdx = 0, pos = 0;
                    
                for (InPath inPath : inSequences.getInPaths()) {
                    
                    if (inPath.getHeader() == "") {
                        
                        pHeader = inPath.getpUId();
                        
                    }else{
                        
                        pHeader = inPath.getHeader();
                        
                    }
                    
                    pathComponents = inPath.getComponents();
                    
                    for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                        
                        uId = std::get<1>(*component);
                        
                        if (std::get<0>(*component) == 'S') {
                        
                            auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                            
                            std::cout<<(*inSegments)[sIdx].getSeqHeader()<<"\t"<<pos;
                            
                            pos += (*inSegments)[sIdx].getInSequence().size();
                            
                            std::cout<<"\t"<<pos<<"\n";
                            
                        }else{
                            
                            auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the segment index
                            
                            pos += (*inGaps)[gIdx].getDist();
                            
                        }
                        
                        
                    }
                    
                    pos = 0;
                    
                }
                
                break;
                
            }
                
            case 'g': { // gaps
                
                std::string pHeader;
                std::vector<InPath> inPaths = inSequences.getInPaths();
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                std::vector<InGap>* inGaps = inSequences.getInGaps();
                std::vector<PathTuple> pathComponents;
                
                unsigned int uId = 0, sIdx = 0, gIdx = 0, pos = 0;
                    
                for (InPath inPath : inSequences.getInPaths()) {
                    
                    if (inPath.getHeader() == "") {
                        
                        pHeader = inPath.getpUId();
                        
                    }else{
                        
                        pHeader = inPath.getHeader();
                        
                    }
                    
                    pathComponents = inPath.getComponents();
                    
                    for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                        
                        uId = std::get<1>(*component);
                        
                        if (std::get<0>(*component) == 'S') {
                        
                            auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                            
                            pos += (*inSegments)[sIdx].getInSequence().size();
                            
                        }else{
                            
                            auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the segment index
                            
                            std::cout<<(*inGaps)[gIdx].getgHeader()<<"\t"<<pos;
                            
                            pos += (*inGaps)[gIdx].getDist();
                            
                            std::cout<<"\t"<<pos<<"\n";
                            
                        }
                        
                    }
                    
                    pos = 0;
                    
                }
                
                break;
                
            }
                
            default:
            case 'a': { // both contigs and gaps in .agp format
                
                std::string pHeader;
                std::vector<InPath> inPaths = inSequences.getInPaths();
                std::vector<InSegment>* inSegments = inSequences.getInSegments();
                std::vector<InGap>* inGaps = inSequences.getInGaps();
                std::vector<PathTuple> pathComponents;
                
                unsigned int uId = 0, sIdx = 0, gIdx = 0, pos = 0;
                    
                for (InPath inPath : inSequences.getInPaths()) {
                    
                    if (inPath.getHeader() == "") {
                        
                        pHeader = inPath.getpUId();
                        
                    }else{
                        
                        pHeader = inPath.getHeader();
                        
                    }
                    
                    pathComponents = inPath.getComponents();
                    
                    for (std::vector<PathTuple>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                        
                        uId = std::get<1>(*component);
                        
                        if (std::get<0>(*component) == 'S') {
                        
                            auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                            
                            std::cout<<pHeader<<"\t"<<pos+1;
                            
                            pos += (*inSegments)[sIdx].getInSequence().size();
                            
                            counter++;
                            
                            std::cout<<"\t"<<pos<<"\t"<<counter<<"\tW\t"<<(*inSegments)[sIdx].getSeqHeader()<<"\t1\t"<<(*inSegments)[sIdx].getInSequence().size()<<"\t"<<std::get<2>(*component)<<"\n";
                            
                        }else{
                            
                            auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                            
                            if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the segment index
                            
                            std::cout<<pHeader<<"\t"<<pos+1;
                            
                            pos += (*inGaps)[gIdx].getDist();
                            
                            counter++;
                            
                            std::cout<<"\t"<<pos<<"\t"<<counter<<"\tN\t"<<(*inGaps)[gIdx].getDist()<<"\t"<<(*inGaps)[gIdx].getgHeader()<<"\tyes\n";
                            
                        }
                        
                    }
                    
                    pos = 0;
                    counter = 0;
                    
                }
            
            }
                
        }
        
        return true;
        
    }
    
    bool reportStats(InSequences &inSequences, unsigned long long int gSize, int bedOutType) { // method to output all summary statistics for the entire sequence set
        std::cout << std::fixed; // disables scientific notation
        std::cout << std::setprecision(2); // 2 decimal poinst

        if (!tabular_flag) {
        
            std::cout<<output("+++Summary+++")<<"\n";
        
        }
        
        if (gSize > 0) {
        
            std::cout<<output("Expected genome size")<<gSize<<"\n";
        
        }
        
        std::cout<<output("# scaffolds")<<inSequences.getScaffN()<<"\n";
        std::cout<<output("Total scaffold length")<<inSequences.getTotScaffLen()<<"\n";
        std::cout<<output("Average scaffold length") << gfa_round(inSequences.computeAverageScaffLen()) << "\n";
        inSequences.evalNstars('s', gSize); // scaffold N* statistics
        std::cout<<output("Scaffold N50")<<inSequences.getScaffN50()<<"\n";
        inSequences.evalAuN('s', gSize); // scaffold auN
        std::cout<<output("Scaffold auN") << gfa_round(inSequences.getScaffauN()) << "\n";
        std::cout<<output("Scaffold L50")<<inSequences.getScaffL50()<<"\n";
        
        if (gSize > 0) {
            
            std::cout<<output("Scaffold NG50")<<inSequences.getScaffNG50()<<"\n";
            std::cout<<output("Scaffold auNG") << gfa_round(inSequences.getScaffauNG()) << "\n";
            std::cout<<output("Scaffold LG50")<<inSequences.getScaffLG50()<<"\n";
            
        }
        std::cout<<output("Largest scaffold")<<inSequences.getLargestScaffold()<<"\n";
        
        std::cout<<output("# contigs")<<inSequences.getSegmentN()<<"\n";
        std::cout<<output("Total contig length")<<inSequences.getTotSegmentLen()<<"\n";
        std::cout<<output("Average contig length") << gfa_round(inSequences.computeAverageSegmentLen()) << "\n";
        inSequences.evalNstars('c', gSize); // contig N* statistics
        std::cout<<output("Contig N50")<<inSequences.getContigN50()<<"\n";
        inSequences.evalAuN('c', gSize); // contig auN
        std::cout<<output("Contig auN") << gfa_round(inSequences.getContigauN()) << "\n";
        std::cout<<output("Contig L50")<<inSequences.getContigL50()<<"\n";
        
        if (gSize > 0) {
            
            std::cout<<output("Contig NG50")<<inSequences.getContigNG50()<<"\n";
            std::cout<<output("Contig auNG") << gfa_round(inSequences.getContigauNG()) << "\n";
            std::cout<<output("Contig LG50")<<inSequences.getContigLG50()<<"\n";
            
        }
        std::cout<<output("Largest contig")<<inSequences.getLargestContig()<<"\n";
        
        std::cout<<output("# gaps")<<inSequences.getGapN()<<"\n";
        std::cout<<output("Total gap length")<<inSequences.getTotGapLen()<<"\n";
        std::cout<<output("Average gap length") << gfa_round(inSequences.computeAverageGapLen()) << "\n";
        inSequences.evalNstars('g'); // gap N* statistics
        std::cout<<output("Gap N50")<<inSequences.getGapN50()<<"\n";
        inSequences.evalAuN('g'); // gap auN
        std::cout<<output("Gap auN") << gfa_round(inSequences.getGapauN()) << "\n";
        std::cout<<output("Gap L50")<<inSequences.getGapL50()<<"\n";
        std::cout<<output("Largest gap")<<inSequences.getLargestGap()<<"\n";
        
        std::cout<<output("Base composition (A:C:G:T)");
        std::cout << inSequences.getTotA() << ":"
                  << inSequences.getTotC() << ":"
                  << inSequences.getTotG() << ":"
                  << inSequences.getTotT() << "\n";
        std::cout<<output("GC content %") << gfa_round(inSequences.computeGCcontent()) << "\n";
        std::cout<<output("# soft-masked bases")<<inSequences.getTotLowerCount()<<"\n";
        
        counter = 0;
        unsigned int connectedComponents = 0;
        
        unsigned int edgeN = inSequences.getEdgeN();

        if (edgeN > 0) {
            
            std::cout<<output("# edges")<<edgeN<<"\n";
            std::cout<<output("Average degree")<<(double)inSequences.getEdgeN()/inSequences.getSegmentN()<<"\n";
        
            inSequences.buildEdgeGraph(inSequences.getEdges());

            verbose("Graph DFS");
            
            std::vector<InSegment>* inSegments = inSequences.getInSegments();
            std::vector<unsigned int> componentLengths;
            unsigned int componentLength = 0;
            
            for (InSegment inSegment : *inSegments) { // loop through all nodes
                
                if (!inSequences.getVisited(inSegment.getuId()) && !inSequences.getDeleted(inSegment.getuId())) { // check if the node was already visited
                    
                    inSequences.dfsEdges(inSegment.getuId(), &componentLength); // if not, visit all connected components recursively
                    connectedComponents++;
                    componentLengths.push_back(componentLength);
                    componentLength = 0;

                }
                
            }

            sort(componentLengths.begin(), componentLengths.end(), std::greater<unsigned int>());

            std::cout<<output("# connected components")<<connectedComponents-inSequences.getDisconnectedComponents()<<"\n";
            std::cout<<output("Largest connected component length")<<componentLengths[0]<<"\n";
            std::cout<<output("# dead ends")<<inSequences.getDeadEnds()<<"\n";
            std::cout<<output("# disconnected components")<<inSequences.getDisconnectedComponents()<<"\n";
            std::cout<<output("Total length disconnected components")<<inSequences.getLengthDisconnectedComponents()<<"\n";
            std::cout<<output("# separated components")<<connectedComponents<<"\n";
        }
        
        unsigned int pathN = inSequences.getPathN();
        
        if (pathN > 0) {
            
            std::cout<<output("# paths")<<pathN<<"\n";
            
        }

        return true;
        
    }
    
    bool nstarReport(InSequences &inSequences, unsigned long long int gSize) { // method to generate all N** reports
        std::cout << std::fixed; // disables scientific notation
        std::cout << std::setprecision(2); // 2 decimal poinst

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

#endif /* GFASTATS_OUTPUT_H */
