//
//  gfastats-input.h
//  gfastats
//
//  Created by Giulio Formenti on 1/16/22.
//

#ifndef GFASTATS_INPUT_H
#define GFASTATS_INPUT_H

class InFile {
    
    std::string h;
    char* c;
    
public:
    
    void readFiles(InSequences &inSequences, std::string &iSeqFileArg, std::string &iSakFileArg, std::string &iAgpFileArg, std::string &iBedIncludeFileArg, std::string &iBedExcludeFileArg, BedCoordinates &bedIncludeList, bool isPipe, char &pipeType, std::string  sortType) {
        
        std::string newLine, seqHeader, seqComment, inSequence, inSequenceQuality, line, bedHeader;
        
        std::unique_ptr<std::istream> stream;
        
        std::vector<Instruction> instructions;
        
        unsigned int begin = 0, end = 0;
 
        if (!iSakFileArg.empty() || (isPipe && (pipeType == 'k'))) {
            
            if (isPipe && (pipeType == 'k')) {
                
                stream = make_unique<std::istream>(std::cin.rdbuf());
                
            }else{
                
                stream = make_unique<std::ifstream>(std::ifstream(iSakFileArg));
                
            }
            
            SAK sak; // create a new swiss army knife
            
            while (getline(*stream, line)) {
                
                std::istringstream iss(line);
                
                instructions.push_back(sak.readInstruction(line)); // use the swiss army knife to read the instruction
                
            }
            
        }
        
        if (!iBedIncludeFileArg.empty() || (isPipe && (pipeType == 'i'))) {
            
            if (isPipe && (pipeType == 'i')) {
                
                stream = make_unique<std::istream>(std::cin.rdbuf());
                
            }else{
                
                stream = make_unique<std::ifstream>(std::ifstream(iBedIncludeFileArg));
                
            }
            
            while (getline(*stream, line)) {
                
                std::istringstream iss(line);
                iss >> bedHeader >> begin >> end;
                
                bedIncludeList.pushCoordinates(bedHeader, begin, end);
                begin = 0, end = 0;
                
            }
            
        }
        
        BedCoordinates bedExcludeList;
        
        if (!iBedExcludeFileArg.empty() || (isPipe && (pipeType == 'e'))) {
            
            if (isPipe && (pipeType == 'e')) {
                
                stream = make_unique<std::istream>(std::cin.rdbuf());
                
            }else{
                
                stream = make_unique<std::ifstream>(std::ifstream(iBedExcludeFileArg));
                
            }
            
            while (getline(*stream, line)) {
                
                std::istringstream iss(line);
                iss >> bedHeader >> begin >> end;
                
                bedExcludeList.pushCoordinates(bedHeader, begin, end);
                begin = 0, end = 0;
                
            }
            
        }
        
        // stream read variable definition
        std::string firstLine;
        unsigned char buffer;
        bool stopStream = false, updateStats = false, isGzip = false;
        unsigned int seqPos = 0; // to keep track of the original sequence order
        
        // start streaming
        stream = make_unique<std::ifstream>(std::ifstream(iSeqFileArg));

        stream->read((char*)(&buffer), 1);

        if (buffer == 0x1f && (stream->peek() == 0x8b)) { // check if pipe is gzipped

            isGzip = true;

        }

        stream->unget();
        
        std::ifstream is(iSeqFileArg);
        
        // this stream takes input from a gzip compressed file
        zstream::igzstream zfin(is);
        
        if (isGzip) {
                
            stream = make_unique<std::istream>(zfin.rdbuf());
            
        
        }
        
        if (isPipe && (pipeType == 's')) { // input is from pipe
            
            std::cin.read((char*)(&buffer), 1);

            if (buffer == 0x1f && (std::cin.peek() == 0x8b)) { // check if pipe is gzipped

                // this stream takes input from a gzip compressed pipe
                zstream::igzstream zin(std::cin);
                
                stream = make_unique<std::istream>(zin.rdbuf());
                
                std::cout<<"Gz pipe currently not supported\n";
                
                exit(1);

            }else {

                stream = make_unique<std::istream>(std::cin.rdbuf());

            }
            
            stream->unget();

        }
        
        if (stream) {
            
            switch (stream->peek()) {
                    
                case '>': {
                    
                    stream->get();
                    
                    while (!stream->eof()) {
                        
                        if(bedIncludeList.size() - bedExcludeList.size() != 0 && bedIncludeList.size() - bedExcludeList.size() == inSequences.getPathN()) { // we have all the sequences needed
                            verbose("Found all sequences, stop streaming input");
                            break;
                        
                        }
                        
                        getline(*stream, newLine);
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }
                        
                        inSequence.clear();
                        
                        getline(*stream, inSequence, '>');
                        
                        inSequence.erase(std::remove(inSequence.begin(), inSequence.end(), '\n'), inSequence.end());
                        
                        verbose("Individual fasta sequence read");
                        
                        Sequence sequence = includeExcludeSeq(seqHeader, seqComment, inSequence, bedIncludeList, bedExcludeList);
                        
                        if (sequence.header != "") {
                            
                            sequence.seqPos = seqPos; // remember the order
                            
                            inSequences.appendSequence(sequence);
                            
                            seqPos++;
                            
                        }
                        
                    }
                    
                    inSequences.joinThreads();
                    
                    break;
                }
                case '@': {
                    
                    while (getline(*stream, newLine)) { // file input
                        
                        newLine.erase(0, 1);
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }
                        
                        getline(*stream, newLine);
                        inSequence = newLine;
                        
                        getline(*stream, newLine);
                        
                        getline(*stream, newLine);
                        inSequenceQuality = newLine;

                        Sequence sequence = includeExcludeSeq(seqHeader, seqComment, inSequence, bedIncludeList, bedExcludeList, inSequenceQuality);
                        
                        if (sequence.header != "") {
                        
                            inSequences.appendSequence(sequence);
                            
                        }
                        
                    }
                    
                    break;
                    
                }
                default: {
                    
                    std::string h_col1, h_col2, h_col3, s, version, gHeader, eHeader, cigar, startS, endS;
                    char sId1Or, sId2Or;
                    
                    InGap gap;
                    InEdge edge;
                    InPath path;
                    unsigned int sId1 = 0, sId2 = 0, dist = 0, start = 0, end = 0;
                    phmap::flat_hash_map<std::string, unsigned int> hash;
                    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got;
                    
                    unsigned int lineN = 1;
                    unsigned int uId = 0, guId = 0, euId = 0;
                    
                    std::string delimiter = "\t";
                    std::vector<std::string> arguments, components; // process the columns of each row
                    
                    size_t pos = 0;
                    
                    getline(*stream, newLine);
                    
                    firstLine = newLine;
                    
                    while ((pos = firstLine.find(delimiter)) != std::string::npos) { // process first line
                        
                        arguments.push_back(firstLine.substr(0, pos));
                        
                        firstLine.erase(0, pos + delimiter.length());
                    
                    }
                    
                    arguments.push_back(firstLine); // last column
                    
                    if (arguments[0] == "H") {
                        
                        h_col2 = arguments[1]; // read header col2
                        delimiter = ":";
                        
                        arguments.clear();
                        
                        while ((pos = h_col2.find(delimiter)) != std::string::npos) { // process version tag
                            
                            arguments.push_back(h_col2.substr(0, pos));
                            
                            h_col2.erase(0, pos + delimiter.length());
                        
                        }
                        
                        arguments.push_back(h_col2); // last column
                        
                        if (arguments[2] != "") {
                            
                            version = arguments[2];
                            verbose("GFA version: " + version);
                            
                        }else{
                            
                            verbose("Failed to parse GFA version. Please check your header.");
                            
                        }
                    
                    }else{
                            
                        verbose("Cannot recognize GFA version from first line. Trying to detect from content.");
                        
                        if (arguments[0] == "S") {
                            
                            if (isInt(arguments[2]) || arguments[2] == "*") {
                                
                                version = '2';
                                verbose("Proposed GFA version: " + version);
                                
                            }else{
                                
                                version = '1';
                                verbose("Proposed GFA version: " + version);
                                
                            }
                            
                        }else if (arguments[0] == "G" || arguments[0] == "O") {
                            
                            version = '2';
                            verbose("Proposed GFA version: " + version);
                            
                        }
                            
                    }
                    
                    if (version[0] == '2') {
                    
                        while (getline(*stream, newLine)) {
                            
                            if (stopStream) {break;}
                            
                            switch (newLine[0]) {
                                    
                                case 'S': {
                                    
                                    strtok(strdup(newLine.c_str()),"\t"); //process first line
                                    seqHeader = strtok(NULL,"\t");
                                    
                                    strtok(NULL,"\t");
                                    s = strtok(NULL,"\t");
                                    inSequence = s;
                                    
                                    c = strtok(NULL,"\t");
                                    if (c != NULL) {
                                        
                                        seqComment = std::string(c);
                                        
                                    }
                                    
                                    stopStream = includeExcludeAppendSegment(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList);
                                    
                                    lineN++;
                                    
                                    break;
                                    
                                }
                                case 'G': {
                                    
                                    delimiter = "\t";

                                    arguments.clear();
                                    
                                    while ((pos = newLine.find(delimiter)) != std::string::npos) {
                                        
                                        arguments.push_back(newLine.substr(0, pos));
                                        
                                        newLine.erase(0, pos + delimiter.length());
                                    
                                    }
                                    
                                    arguments.push_back(newLine); // last column
                                    
                                    gHeader = arguments[1];
                                    
                                    uId = inSequences.getuId();
                                    
                                    inSequences.insertHash(gHeader, uId);
                                    
                                    guId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                    
                                    sId1Or = arguments[2].back(); // get sequence orientation in the gap
                                    
                                    seqHeader = std::string(arguments[2]);
                                    seqHeader.pop_back();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId1 = uId;
                                        
                                        inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                        
                                    }else{
                                        
                                        sId1 = got->second;
                                        
                                    }
                                    
                                    sId2Or = arguments[3].back(); // get sequence orientation in the gap
                                    
                                    seqHeader = arguments[3];
                                    seqHeader.pop_back();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId2 = uId;
                                        
                                        inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                        
                                    }else{
                                        
                                        sId2 = got->second;
                                        
                                    }
                                    
                                    dist = stoi(arguments[4]);
                                    
                                    verbose("Processing gap " + gHeader + " (uId: " + std::to_string(uId) + ")");
                                    
                                    gap.newGap(guId, sId1, sId2, sId1Or, sId2Or, dist, gHeader);
                                    
                                    inSequences.addGap(gap);
                                    
                                    lineN++;
                                                 
                                    break;
                                    
                                }

                                case 'E': {
                             
                                    delimiter = "\t";

                                    arguments.clear();
                                    
                                    while ((pos = newLine.find(delimiter)) != std::string::npos) {
                                        
                                        arguments.push_back(newLine.substr(0, pos));
                                        
                                        newLine.erase(0, pos + delimiter.length());
                                    
                                    }
                                    
                                    arguments.push_back(newLine); // last column
                                    
                                    eHeader = arguments[1];
                                    
                                    uId = inSequences.getuId();
                                    
                                    inSequences.insertHash(eHeader, uId);
                                    
                                    euId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                    
                                    sId1Or = arguments[2].back(); // get sequence orientation in the edge
                                    
                                    seqHeader = std::string(arguments[2]);
                                    seqHeader.pop_back();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId1 = uId;
                                        
                                        inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                        
                                    }else{
                                        
                                        sId1 = got->second;
                                        
                                    }
                                    
                                    sId2Or = arguments[3].back(); // get sequence orientation in the edge
                                    
                                    seqHeader = arguments[3];
                                    seqHeader.pop_back();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId2 = uId;
                                        
                                        inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                        
                                    }else{
                                        
                                        sId2 = got->second;
                                        
                                    }                            
      
                                    cigar = arguments[8];  
                                    
                                    edge.newEdge(euId, sId1, sId2, sId1Or, sId2Or, cigar, eHeader);
                                    
                                    inSequences.appendEdge(edge);
 
                                    lineN++;
                                                 
                                    break;
                                    
                                }
                                    
                                case 'O': {
                                    
                                    arguments.clear();
                                    
                                    arguments = readDelimited(newLine, "\t");
                                    
                                    seqHeader = arguments[1];
                                    
                                    uId = inSequences.getuId();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table to look for the header
                                    
                                    if (got == hash.end()) { // this is the first time we see this header
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                        
                                    }else{
                                        
                                        fprintf(stderr, "Error: path name already exists (%s). Terminating.\n", seqHeader.c_str()); exit(1);
                                        
                                    }
                                    
                                    path.newPath(uId, seqHeader);
                                    
                                    inSequences.uId.next();
                                    
                                    components = readDelimited(arguments[2], " ");
                                    
                                    for (std::string component : components) {
                                        
                                        sId1Or = component.back(); // get sequence orientation
                                        
                                        if (sId1Or == '+' || sId1Or == '-') { // only segments have orientation
                                        
                                            component.pop_back();
                                            
                                        }
                                        
                                        if (component.find("(") != std::string::npos && component.find(":") != std::string::npos && component.find(")") != std::string::npos) {
                                            
                                            startS = component.substr(component.find("(") + 1, component.find(":") - component.find("(") - 1);
                                            endS = component.substr(component.find(":") + 1, component.find(")") - component.find(":") - 1);
                                            
                                            start = std::stoi(startS);
                                            end = std::stoi(endS);

                                            component = component.substr(0, component.find("("));
                                            
                                        }else{
                                            
                                            start = 0;
                                            end = 0;
                                            
                                        }
                                        
                                        if (end != 0) {
                                        
                                            verbose("Adding only coordinates " + std::to_string(start) + ":" + std::to_string(end) + "(" + component + ")");
                                            
                                        }
                                    
                                        hash = inSequences.getHash1();
                                        
                                        got = hash.find(component); // get the headers to uIds table (remove sequence orientation in the gap first)
                                        
                                        if (got == hash.end()) { // this is the first time we see this segment
                                            
                                            uId = inSequences.getuId();
                                            
                                            inSequences.insertHash(component, uId);
                                        
                                            sId1 = uId;
                                            
                                            inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                            
                                        }else{
                                            
                                            sId1 = got->second;
                                            
                                        }
                                        
                                        std::vector<InSegment>* inSegments = inSequences.getInSegments();
                                        std::vector<InGap>* inGaps = inSequences.getInGaps();
                                        
                                        auto sId = find_if(inSegments->begin(), inSegments->end(), [sId1](InSegment& obj) {return obj.getuId() == sId1;}); // given a uId, find it in nodes
                                    
                                        if (sId != inSegments->end()) {
                                            
                                            path.add(SEGMENT, sId1, sId1Or, start, end);
                                             
                                        }else{
                                            
                                            auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1](InGap& obj) {return obj.getuId() == sId1;}); // given a uId, find it in gaps
                                            
                                            if (gId != inGaps->end()) {
                                            
                                                path.add(GAP, sId1, '0', start, end);
                                            
                                            }
                                            
                                        }
                                        
                                    }
                                    
                                    c = strtok(NULL,"\t");
                                    if (c != NULL) {
                                        
                                        seqComment = std::string(c);
                                        path.setComment(seqComment);
                                        
                                    }
                                    
                                    inSequences.addPath(path);
                                    
                                    lineN++;
                                    
                                    break;
                                    
                                }
                                    
                            }
                            
                        }
                        
                    }else if (version[0] == '1') {
                    
                        while (getline(*stream, newLine)) {
                            
                            if (stopStream) {break;}
                            
                            switch (newLine[0]) {
                                    
                                case 'S': {
                                    
                                    strtok(strdup(newLine.c_str()),"\t"); // process first line
                                    h = strtok(NULL,"\t");
                                    
                                    seqHeader = h;
                                    
                                    s = strtok(NULL,"\t");
                                    inSequence = s;
                                    
                                    c = strtok(NULL,"\t");
                                    if (c != NULL) {
                                        
                                        seqComment = std::string(c);
                                        
                                    }
                                    
                                    stopStream = includeExcludeAppendSegment(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList);
                                    
                                    lineN++;
                                    
                                    break;
                                    
                                }

                                case 'L': {
                                    
                                    delimiter = "\t";
                                    
                                    arguments.clear();
                                    
                                    while ((pos = newLine.find(delimiter)) != std::string::npos) {
                                        
                                        arguments.push_back(newLine.substr(0, pos));
                                        
                                        newLine.erase(0, pos + delimiter.length());
                                    
                                    }

                                    uId = inSequences.getuId();
                                    
                                    euId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                    
                                    arguments.push_back(newLine); // last column
                                    
                                    sId1Or = arguments[2][0]; // get sequence orientation in the edge
                                    
                                    seqHeader = arguments[1];
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId1 = uId;
                                        
                                        uId++;
                                        
                                        inSequences.uId.next(); // we have touched a segment need to increase the unique segment counter
                                        
                                    }else{
                                        
                                        sId1 = got->second;
                                        
                                    }
                                    
                                    sId2Or = arguments[4][0]; // get sequence orientation in the edge
                                    
                                    seqHeader = arguments[3];
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId2 = uId;
                                        
                                        uId++;
                                        
                                        inSequences.uId.next(); // we have touched a segment need to increase the unique segment counter
                                        
                                    }else{
                                        
                                        sId2 = got->second;
                                        
                                    }
                                    
                                    cigar = arguments[5];  
                                    
                                    edge.newEdge(euId, sId1, sId2, sId1Or, sId2Or, cigar);
                                    
                                    inSequences.appendEdge(edge);
 
                                    lineN++;
                                                 
                                    break;
                                    
                                }

                                case 'P': {
                                    
                                    strtok(strdup(newLine.c_str()),"\t"); // process first line
                                    
                                    seqHeader = strtok(NULL,"\t");

                                    uId = inSequences.getuId();
                                    inSequences.uId.next();
                                    
                                    path.newPath(uId, seqHeader);
                                    
                                    s = strtok(NULL,"\t");
                                    
                                    delimiter = ",";
                                    
                                    arguments.clear();
                                    
                                    while ((pos = s.find(delimiter)) != std::string::npos) {
                                        
                                        arguments.push_back(s.substr(0, pos));
                                        
                                        s.erase(0, pos + delimiter.length());
                                    
                                    }
                                    
                                    arguments.push_back(s); // last column
                                    
                                    for (std::string component : arguments) {
                                        
                                        sId1Or = component.back(); // get sequence orientation
                                        
                                        if (sId1Or == '+' || sId1Or == '-') { // only segments have orientation
                                        
                                            component.pop_back();
                                            
                                        }
                                    
                                        hash = inSequences.getHash1();
                                        
                                        got = hash.find(component); // get the headers to uIds table (remove sequence orientation in the gap first)
                                        
                                        if (got == hash.end()) { // this is the first time we see this segment
                                            
                                            uId = inSequences.getuId();
                                            
                                            inSequences.insertHash(component, uId);
                                        
                                            sId1 = uId;
                                            
                                            inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                            
                                        }else{
                                            
                                            sId1 = got->second;
                                            
                                        }
                                        
                                        std::vector<InSegment>* inSegments = inSequences.getInSegments();
                                        std::vector<InGap>* inGaps = inSequences.getInGaps();
                                        
                                        auto sId = find_if(inSegments->begin(), inSegments->end(), [sId1](InSegment& obj) {return obj.getuId() == sId1;}); // given a uId, find it in nodes
                                    
                                        if (sId != inSegments->end()) {
                                            
                                            path.add(SEGMENT, sId1, sId1Or);
                                             
                                        }else{
                                            
                                            auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1](InGap& obj) {return obj.getuId() == sId1;}); // given a uId, find it in gaps
                                            
                                            if (gId != inGaps->end()) {
                                            
                                                path.add(GAP, sId1, '0');
                                            
                                            }
                                            
                                        }
                                        
                                    }
                                    strtok(NULL,"\t");
                                    
                                    c = strtok(NULL,"\t");
                                    
                                    if (c != NULL) {
                                        
                                        seqComment = std::string(c);
                                        path.setComment(seqComment);
                                        
                                    }
                                    
                                    inSequences.addPath(path);
                                    
                                    lineN++;
                                    
                                    break;
                                    
                                }
                                    
                            }
                            
                        }
                        
                        break;
                        
                    }
                    
                }
                
            }
            
            verbose("End of file");
                
        }else{

            fprintf(stderr, "Stream not successful: %s", iSeqFileArg.c_str());
            exit(1);

        }
        
        if (rmGaps_flag) {
         
            inSequences.removeTerminalGaps();
            
        }
        
        if (discoverPaths_flag) {
            
            inSequences.discoverPaths();
            
        }
        
        if (!instructions.empty()) {
            
            verbose("\nStarted instruction execution");
        
            SAK sak; // create a new swiss army knife
            
            for (Instruction instruction : instructions) { // execute swiss army knife instructions
                
                sak.executeInstruction(inSequences, instruction);
                
                verbose(instruction.action + " instruction executed");
                
            }
        
        }
        
        inSequences.updateStats();

        if (!iAgpFileArg.empty() || (isPipe && (pipeType == 'a'))) {
            
            std::string pHeaderNew, pHeader1, pHeader2, gHeader, instruction, coord1, coord2;
            char pId1Or = '+', pId2Or;
            
            unsigned int pUId = 0, pUId1 = 0, pUId2 = 0, gUId = 0, dist = 0, seqLen, pathLen, start1 = 0, end1 = 0, start2 = 0, end2 = 0;
            phmap::flat_hash_map<std::string, unsigned int> hash;
            phmap::flat_hash_map<std::string, unsigned int>::const_iterator got;
            
            std::vector<std::string> arguments; // line arguments
            std::vector<unsigned int> oldPaths; // vector of paths flagged to be removed only at the end
            
            std::vector<InPath> paths = inSequences.getInPaths();
            
            for (InPath path : paths) {
                
                oldPaths.push_back(path.getpUId());
                
            }
            
            if (isPipe && (pipeType == 'a')) {
                
                stream = make_unique<std::istream>(std::cin.rdbuf());
                
            }else{
                
                stream = make_unique<std::ifstream>(std::ifstream(iAgpFileArg));
                
            }

            std::queue<std::string> nextLines;
            
            while (true) {
                if(nextLines.size() > 0) {
                    line = nextLines.front();
                    nextLines.pop();
                } else if(!getline(*stream, line)) {
                    break;
                }
                
                std::istringstream iss(line); // line to string
                
                arguments = readDelimited(line, "\t", "#"); // read the columns in the line
                
                if (arguments.size() == 0) {continue;}
                
                pHeaderNew = arguments[0]; // this is the current header
                
                if (arguments[4] == "W") { // this is an old path
                    
                    if (!discoverPaths_flag) {
                    
                        pHeader1 = arguments[5];
                    
                    }else{
                        
                        pHeader1 = arguments[5] + "_path";
                    }
                    
                    pId1Or = arguments[8][0];
                    
                    hash = inSequences.getHash1();
                    
                    got = hash.find(pHeader1); // get the headers to uIds table
                    
                    if (got != hash.end()) { // this is not the first time we see this path
                        
                        pUId1 = got->second;
                        
                    }else{
                        
                        fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader1.c_str()); // sequence not found
                        
                        continue;
                        
                    }
                    
                    pathLen = inSequences.pathLen(pUId1);
                    
                    start1 = stoi(arguments[6]);
                    end1 = stoi(arguments[7]);
                    
                    seqLen = end1 - start1 + 1;
                    
                    if(seqLen != pathLen) {

                        fprintf(stderr, "Warning: sequence length (%u) differs from path length (%u). Subsetting (%s).\n", seqLen, pathLen, pHeader1.c_str());

                    }else{
                        
                        start1 = 0;
                        end1 = 0;
                        
                    }
                    
                    getline(*stream, line);
                    nextLines.push(line);
                    
                    arguments = readDelimited(line, "\t", "#"); // read the next sequence
                    
                    if(pHeaderNew != arguments[0]) { // if this path does not need to be joined to anything that follows, we create a new path
                        
                        InPath path;
                        
                        pUId = inSequences.getuId();
                        
                        path.newPath(pUId, pHeaderNew);
                        
                        std::vector<PathComponent> pathComponents = inSequences.getInPath(pUId1).getComponents();
                        
                        path.append({std::begin(pathComponents), std::end(pathComponents)});
                        
                        inSequences.insertHash(pHeaderNew, pUId);
                        
                        inSequences.uId.next();
                        
                        inSequences.addPath(path);
                        
                        if(seqLen != pathLen) { // if it also needs to be trimmed
                            
                            inSequences.trimPathByUId(pUId, start1, end1);
                            
                        }
                        
                        if(pId1Or == '-') {
                            
                            inSequences.revComPath(pUId);
                            
                        }
                        
                    }
                
                    
                }else if(arguments[4] == "N" || arguments[4] == "U"){

                    hash = inSequences.getHash1();
                    
                    got = hash.find(pHeader1); // get the headers to uIds table (remove sequence orientation in the gap first)
                    
                    if (got == hash.end()) { // this is the first time we see this path
                        
                        fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader1.c_str()); // if the preceding sequence was not found we do not introduce a gap
                        
                        continue;
                        
                    }
                    
                    gUId = inSequences.getuId();
                    
                    if (arguments[6] == "scaffold") {
                        
                        hash = inSequences.getHash1();
                        
                        got = hash.find("gap"+std::to_string(gUId)); // get the headers to uIds table
                        
                        while (got != hash.end()) { // this is not the first time we see this path
                            
                            gUId++;
                            
                            got = hash.find("gap"+std::to_string(gUId)); // get the headers to uIds table
                            
                            inSequences.uId.next();
                        }
                    
                        inSequences.uId.next();
                        gHeader = "gap"+std::to_string(gUId);
                    
                    }else{
                        
                        gHeader = arguments[6];
                        inSequences.uId.next();
                        
                    }
                    
                    inSequences.insertHash(gHeader, gUId);
                    
                    
                    dist = stoi(arguments[5]);
                    
                    getline(*stream, line);
                    
                    arguments = readDelimited(line, "\t", "#"); // read the next sequence
                    
                    if (arguments.size() == 0) {continue;}
                    
                    if (!discoverPaths_flag) {
                    
                        pHeader2 = arguments[5];
                    
                    }else{
                        
                        pHeader2 = arguments[5] + "_path";
                    }
                        
                    pId2Or = arguments[8][0];
                    
                    hash = inSequences.getHash1();
                    
                    got = hash.find(pHeader2); // get the headers to uIds table (remove sequence orientation in the gap first)
                    
                    if (got != hash.end()) { // this is not the first time we see this path
                        
                        pUId2 = got->second;
                        
                    }else{
                        
                        fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader2.c_str()); // sequence not found
                        
                        continue;
                        
                    }
                    
                    pathLen = inSequences.pathLen(pUId2);
                    
                    start2 = stoi(arguments[6]);
                    end2 = stoi(arguments[7]);
                    
                    seqLen = end2 - start2 + 1;
                    
                    if(seqLen != pathLen) {

                        fprintf(stderr, "Warning: sequence length (%u) differs from path length (%u). Subsetting (%s).\n", seqLen, pathLen, pHeader2.c_str());

                    }else{
                        
                        start2 = 0;
                        end2 = 0;
                        
                    }
                    
                    SAK sak;
                    
                    coord1 = start1 != 0 ? "(" + std::to_string(start1) + ":" + std::to_string(end1) + ")" : "";
                    coord2 = start2 != 0 ? "(" + std::to_string(start2) + ":" + std::to_string(end2) + ")" : "";
                    
                    instruction = "JOIN\t" + pHeader1 + coord1 + pId1Or + "\t" + pHeader2 + coord2 + pId2Or + "\t" + std::to_string(dist) + "\t" + gHeader + "\t" + pHeaderNew + "\t" + std::to_string(gUId);
                    
                    fprintf(stderr, "%s\n", instruction.c_str());
                    
                    sak.executeInstruction(inSequences, sak.readInstruction(instruction));
                    
                    pHeader1 = pHeaderNew;
                    start1 = 0;
                    end1 = 0;
                    pId1Or = '+';

                }
                
            }
                
            for (unsigned int pUId : oldPaths) { // remove paths left
                
                inSequences.removePath(pUId, true); // silently remove the original paths that were not joined or duplicated
                
            }
            
            updateStats = true;
            
        }
            
        if (sortType == "ascending") {
            
            inSequences.sortPathsByNameAscending();
            
        }else if (sortType == "descending") {
            
            inSequences.sortPathsByNameDescending();
            
        }else if (sortType == "largest") {
            
            inSequences.sortPathsBySize(0);

        }else if (sortType == "smallest") {
            
            inSequences.sortPathsBySize(1);
            
        }else if (sortType != "none" && ifFileExists(sortType.c_str())){
                
            stream = make_unique<std::ifstream>(std::ifstream(sortType));
            
            std::string header;
            std::vector<std::string> headerList;
            
            while (getline(*stream, line)) { // read the file to vector
                
                std::istringstream iss(line);
                iss >> header;
                
                headerList.push_back(header);
                
            }
            
            inSequences.sortPathsByList(headerList);
            
        }else{
            
            inSequences.sortPathsByOriginal();
            
            
        }
        
        if (updateStats) {
            
            inSequences.updateStats();
            
        }
    }
    
    Sequence includeExcludeSeq(std::string seqHeader, std::string seqComment, std::string inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string inSequenceQuality = "") {
        
        std::vector<std::string> bedIncludeListHeaders;
        std::vector<std::string> bedExcludeListHeaders;
        unsigned int pos = 0, cBegin = 0, cEnd = 0, offset = 0, prevCEnd = 0;

        bedIncludeListHeaders = bedIncludeList.getSeqHeaders();
        bedExcludeListHeaders = bedExcludeList.getSeqHeaders();
        bool outSeq = false;
        
        verbose("Processing sequence: " + seqHeader);
        
        if   (bedIncludeList.empty() &&
              bedExcludeList.empty()) {
            
            outSeq = true;
            
        }else if(!bedIncludeList.empty() &&
                  bedExcludeList.empty()) {
            
            offset = 0, prevCEnd = 0;
            outSeq = false;
            
            auto it = begin(bedIncludeListHeaders);
            
            while (it != end(bedIncludeListHeaders)) {
                
                it = std::find(it, bedIncludeListHeaders.end(), seqHeader);
                
                if (it == end(bedIncludeListHeaders) || bedIncludeList.getSeqHeader(pos) != seqHeader) {
                    
                    break;
                    
                }
                
                outSeq = true;

                cBegin = bedIncludeList.getcBegin(pos);
                cEnd = bedIncludeList.getcEnd(pos);
                
                if (!(cBegin == 0 && cEnd == 0)) {
                    
                    inSequence.erase(offset, cBegin-prevCEnd);
                    
                    if (inSequenceQuality != "") {
                    
                        inSequenceQuality.erase(offset, cBegin-prevCEnd);
                    
                    }
                        
                    offset += cEnd-cBegin;
                    prevCEnd = cEnd;
                    
                }
              
                ++it;
                pos++;
                
            }
                
            if (outSeq && inSequence.size()>0) {
                
                if (offset>0) {
                
                    inSequence.erase(offset, inSequence.size()-offset);
                    
                    if (inSequenceQuality != "") {
                    
                        inSequenceQuality.erase(offset, inSequenceQuality.size()-offset);
                        
                    }
                    
                }
                
                outSeq = true;
            
            }
                
        }else if(bedIncludeList.empty() &&
                !bedExcludeList.empty()) {
            
            pos = 0;
            offset = 0;
            outSeq = true;
            
            auto it = begin(bedExcludeListHeaders);
            
            while (it != end(bedExcludeListHeaders)) {
                
                it = std::find(it, bedExcludeListHeaders.end(), seqHeader);
                
                if (it == end(bedExcludeListHeaders)) {
                    
                    break;
                    
                }

                cBegin = bedExcludeList.getcBegin(pos);
                cEnd = bedExcludeList.getcEnd(pos);
                
                if (!(cBegin == 0 && cEnd == 0)) {
                    
                    inSequence.erase(cBegin-offset, cEnd-cBegin);
                    
                    if (inSequenceQuality != "") {
                    
                        inSequenceQuality.erase(cBegin-offset, cEnd-cBegin);
                        
                    }
                        
                    offset += cEnd-cBegin;
                    
                }else{
                    
                    outSeq = false;
                    
                }
              
                ++it;
                pos++;
                
            }
                    
        }else if
                (!bedIncludeList.empty() &&
                 !bedExcludeList.empty() &&
                 std::find(bedIncludeListHeaders.begin(), bedIncludeListHeaders.end(), seqHeader) != bedIncludeListHeaders.end() &&
                 std::find(bedExcludeListHeaders.begin(), bedExcludeListHeaders.end(), seqHeader) == bedExcludeListHeaders.end()) {
                    
                    outSeq = true;
                    
        }
        
        if (outSeq && inSequence.size()>0) {
        
            Sequence sequence {seqHeader, seqComment, inSequence, inSequenceQuality};
            return sequence;
        
        }else {
            
            verbose("Sequence entirely removed as a result of BED filter: " + seqHeader);
            
            Sequence sequence {""};
            return sequence;
            
        }
        
    }
    
    bool includeExcludeAppendSegment(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality = NULL) {
        
        std::vector<std::string> bedIncludeListHeaders;
        std::vector<std::string> bedExcludeListHeaders;
        unsigned int pos = 0, cBegin = 0, cEnd = 0, offset = 0, prevCEnd = 0;
 
        bedIncludeListHeaders = bedIncludeList.getSeqHeaders();
        bedExcludeListHeaders = bedExcludeList.getSeqHeaders();
        bool outSeq;
        
        if   (bedIncludeList.empty() &&
              bedExcludeList.empty()) {
            
            inSequences->appendSegment(seqHeader, seqComment, inSequence, inSequenceQuality);
            
        }else if(!bedIncludeList.empty() &&
                  bedExcludeList.empty()) {
            
            if(inSequences->getInSegments()->size() == bedIncludeList.size()) { // check if we retrieved all we needed
                
                verbose("Found all sequences, stop streaming input");
                
                return true;
                
            }
            
            offset = 0, prevCEnd = 0;
            outSeq = false;
            
            auto it = begin(bedIncludeListHeaders);
            
            while (it != end(bedIncludeListHeaders)) {
                
                it = std::find(it, bedIncludeListHeaders.end(), *seqHeader);
                
                if (it == end(bedIncludeListHeaders) || bedIncludeList.getSeqHeader(pos) != *seqHeader) {
                    
                    break;
                    
                }
                
                outSeq = true;

                cBegin = bedIncludeList.getcBegin(pos);
                cEnd = bedIncludeList.getcEnd(pos);
                
                if (!(cBegin == 0 && cEnd == 0)) {
                    
                    inSequence->erase(offset, cBegin-prevCEnd);
                    
                    if (inSequenceQuality != NULL) {
                    
                        inSequenceQuality->erase(offset, cBegin-prevCEnd);
                    
                    }
                        
                    offset += cEnd-cBegin;
                    prevCEnd = cEnd;
                    
                }
              
                ++it;
                pos++;
                
            }
                
            if (outSeq && inSequence->size()>0) {
                
                if (offset>0) {
                
                    inSequence->erase(offset, inSequence->size()-offset);
                    
                    if (inSequenceQuality != NULL) {
                    
                        inSequenceQuality->erase(offset, inSequenceQuality->size()-offset);
                        
                    }
                    
                }
                
                inSequences->appendSegment(seqHeader, seqComment, inSequence, inSequenceQuality);
            
            }else {
                
                verbose("Scaffold entirely removed as a result of include: " + *seqHeader);
                
            }
                
        }else if(bedIncludeList.empty() &&
                !bedExcludeList.empty()) {
                
            offset = 0;
            outSeq = true;
            
            auto it = begin(bedExcludeListHeaders);
            
            while (it != end(bedExcludeListHeaders)) {
                
                it = std::find(it, bedExcludeListHeaders.end(), *seqHeader);
                
                if (it == end(bedExcludeListHeaders)) {
                    
                    break;
                    
                }

                cBegin = bedExcludeList.getcBegin(pos);
                cEnd = bedExcludeList.getcEnd(pos);
                
                if (!(cBegin == 0 && cEnd == 0)) {
                    
                    inSequence->erase(cBegin-offset, cEnd-cBegin);
                    
                    if (inSequenceQuality != NULL) {
                    
                        inSequenceQuality->erase(cBegin-offset, cEnd-cBegin);
                        
                    }
                        
                    offset += cEnd-cBegin;
                    
                }else{
                    
                    outSeq = false;
                    
                }
              
                ++it;
                pos++;
                
            }
                
            if (outSeq && inSequence->size()>0) {
            
                inSequences->appendSegment(seqHeader, seqComment, inSequence, inSequenceQuality);
            
            }else {
                
                verbose("Scaffold entirely removed as a result of exclude: " + *seqHeader);
                
            }
                    
        }else if
                (!bedIncludeList.empty() &&
                 !bedExcludeList.empty() &&
                 std::find(bedIncludeListHeaders.begin(), bedIncludeListHeaders.end(), *seqHeader) != bedIncludeListHeaders.end() &&
                 std::find(bedExcludeListHeaders.begin(), bedExcludeListHeaders.end(), *seqHeader) == bedExcludeListHeaders.end()) {
                    
                    if(inSequences->getInSegments()->size() == bedIncludeList.size()) { // check if we retrieved all we needed
                        
                        verbose("Found all sequences, stop streaming input");
                        
                        return true;
                        
                    }
                    
                    inSequences->appendSegment(seqHeader, seqComment, inSequence, inSequenceQuality);
                    
        }
        
        return false;
        
    }
    
};

#endif /* GFASTATS_INPUT_H */
