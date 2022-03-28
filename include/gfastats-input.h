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
    std::vector<std::string> bedIncludeListHeaders;
    std::vector<std::string> bedExcludeListHeaders;
    unsigned int pos = 0, cBegin = 0, cEnd = 0, offset = 0, prevCEnd = 0;
    
public:
    
    bool getFasta(std::istream& is, std::string& s) // buffer reader for fasta
    {
        s.clear();
        
//        int i = 0;
//        
//        size_t buffer_s = 500000000; // large input buffer
        
        getline(is, s, '>');
        
        s.erase(std::remove(s.begin(), s.end(), '\n'), s.end());
        
//        char* str_new = (char*) malloc(buffer_s * sizeof(char));
//
//        for(char* c = str; *c != '\0'; c++) { // for each character
//
//            if (*c != '\n') { // remove newline characters
//
//                str_new[i] = *c; i++;
//
//            }
//
//        }

//        s = str_new;
//        free(str);
//        free(str_new);

        return is.eof() ? false : true;

    }
    
    InSequences readFiles(std::string &iSeqFileArg, std::string &iSakFileArg, std::string &iAgpFileArg, std::string &iBedIncludeFileArg, std::string &iBedExcludeFileArg, BedCoordinates &bedIncludeList, bool isPipe, char &pipeType, std::string  sortType) {
        
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
        
        InSequences inSequences;
        
        std::string firstLine;
        char firstChar;
        unsigned char buffer[2];
        bool stopStream;
        
        std::ifstream is(iSeqFileArg);
        
        // this stream takes input from a gzip compressed file
        zstream::igzstream zfin(is);
        
        if (determineGzip(iSeqFileArg)) { // input is a compressed file
            
            stream = make_unique<std::istream>(zfin.rdbuf());
            
        } else if (isPipe && (pipeType == 's')) { // input is from pipe
            
            std::cin.read((char*)(&buffer[0]), 2);

            if (buffer[0] == 0x1f && (buffer[1] == 0x8b)) { // check if pipe is gzipped

                // this stream takes input from a gzip compressed pipe
                zstream::igzstream zin(std::cin);
                
                stream = make_unique<std::istream>(zin.rdbuf());
                
                getline(*stream, newLine);
                
                std::cout<<"Gz pipe currently not supported\n";
                
                exit(1);

            }else{

                stream = make_unique<std::istream>(std::cin.rdbuf());

            }

        } else { // input is a regular plain text file

            stream = make_unique<std::ifstream>(std::ifstream(iSeqFileArg));

        }
        
        if (stream) {
            
            getline(*stream, newLine);
            
            if (isPipe && (pipeType == 's')) { // if input from pipe, recover the first two bytes
                
                newLine.insert (0, 1, buffer[1]);
                newLine.insert (0, 1, buffer[0]);
                
            }
            
            firstLine = newLine;
            firstChar = newLine[0];
            
            switch (firstChar) {
                    
                case '>': {
                        
                    firstLine.erase(0, 1);
                    
                    h = std::string(strtok(strdup(firstLine.c_str())," ")); //process header line
                    c = strtok(NULL,""); //read comment
                    
                    seqHeader = h;
                    
                    if (c != NULL) {
                        
                        seqComment = std::string(c);
                        
                    }
                    
                    while (getFasta(*stream, inSequence)) {
                        
                        verbose("Individual fasta sequence read");
                        
                        stopStream = includeExcludeAppend(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList);
                        
                        getline(*stream, newLine);
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }
                        
                        if (stopStream) {break;}
                        
                    }
                    
                    includeExcludeAppend(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList);
                    
                    break;
                }
                case '@': {
                        
                    firstLine.erase(0, 1);
                    
                    h = std::string(strtok(strdup(firstLine.c_str())," ")); //process header line
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
                    
                    includeExcludeAppend(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList, &inSequenceQuality);
                    
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

                        stopStream = includeExcludeAppend(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList, &inSequenceQuality);
                        
                        if (stopStream) {break;}
                        
                    }
                    
                    break;
                    
                }
                default: {
                    
                    if ((!isPipe || pipeType != 's') && !determineGzip(iSeqFileArg)) {

                        stream->clear();
                        stream->seekg(0, stream->beg);

                    }
                    
                    std::string h_col1, h_col2, h_col3, s, version, gHeader, eHeader, cigar;
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
                            
                        }else if (arguments[0] == "G") {
                            
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
                                    
                                    inSequences.insertHash1(gHeader, uId); // header to hash table
                                    inSequences.insertHash2(uId, gHeader); // uID to hash table
                                    
                                    guId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                    
                                    sId1Or = arguments[2].back(); // get sequence orientation in the gap
                                    
                                    seqHeader = std::string(arguments[2]);
                                    seqHeader.pop_back();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash1(seqHeader, uId); // header to hash table
                                        inSequences.insertHash2(uId, seqHeader); // header to hash table
                                    
                                        sId1 = uId;
                                        
                                        inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                        
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
                                        
                                        inSequences.insertHash1(seqHeader, uId); // header to hash table
                                        inSequences.insertHash2(uId, seqHeader); // header to hash table
                                    
                                        sId2 = uId;
                                        
                                        inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                        
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
                                    
                                    inSequences.insertHash1(eHeader, uId); // header to hash table
                                    inSequences.insertHash2(uId, eHeader); // uID to hash table
                                    
                                    euId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                    
                                    sId1Or = arguments[2].back(); // get sequence orientation in the edge
                                    
                                    seqHeader = std::string(arguments[2]);
                                    seqHeader.pop_back();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash1(seqHeader, uId); // header to hash table
                                        inSequences.insertHash2(uId, seqHeader); // header to hash table
                                    
                                        sId1 = uId;
                                        
                                        inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                        
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
                                        
                                        inSequences.insertHash1(seqHeader, uId); // header to hash table
                                        inSequences.insertHash2(uId, seqHeader); // header to hash table
                                    
                                        sId2 = uId;
                                        
                                        inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                        
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
                                        
                                        inSequences.insertHash1(seqHeader, uId); // header to hash table
                                        inSequences.insertHash2(uId, seqHeader); // header to hash table
                                        
                                    }else{
                                        
                                        fprintf(stderr, "Error: path name already exists (%s). Terminating.\n", seqHeader.c_str()); exit(1);
                                        
                                    }
                                    
                                    path.newPath(uId, seqHeader);
                                    
                                    inSequences.setuId(uId+1);
                                    
                                    components = readDelimited(arguments[2], " ");
                                    
                                    for (std::string component : components) {
                                        
                                        sId1Or = component.back(); // get sequence orientation
                                        
                                        if (sId1Or == '+' || sId1Or == '-') { // only segments have orientation
                                        
                                            component.pop_back();
                                            
                                        }
                                        
                                        if (component.find("(") != std::string::npos && component.find(":") != std::string::npos && component.find(")") != std::string::npos) {
                                            
                                            start = stoi(component.substr(component.find("(") + 1, component.find(":") - component.find("(") - 1));
                                            end = stoi(component.substr(component.find(":") + 1, component.find(")") - component.find(":") - 1));
                                            
                                            component = component.substr(0, component.find("("));
                                            
                                        }else{
                                            
                                            start = 0;
                                            end = 0;
                                            
                                        }
                                        
                                        verbose(component + " " + std::to_string(start) + " " + std::to_string(end));
                                    
                                        hash = inSequences.getHash1();
                                        
                                        got = hash.find(component); // get the headers to uIds table (remove sequence orientation in the gap first)
                                        
                                        if (got == hash.end()) { // this is the first time we see this segment
                                            
                                            uId = inSequences.getuId();
                                            
                                            inSequences.insertHash1(component, uId); // header to hash table
                                            inSequences.insertHash2(uId, component); // uId to hash table
                                        
                                            sId1 = uId;
                                            
                                            inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                            
                                        }else{
                                            
                                            sId1 = got->second;
                                            
                                        }
                                        
                                        std::vector<InSegment>* inSegments = inSequences.getInSegments();
                                        std::vector<InGap>* inGaps = inSequences.getInGaps();
                                        
                                        auto sId = find_if(inSegments->begin(), inSegments->end(), [sId1](InSegment& obj) {return obj.getuId() == sId1;}); // given a uId, find it in nodes
                                    
                                        if (sId != inSegments->end()) {
                                            
                                            path.add('S', sId1, sId1Or, start, end);
                                             
                                        }else{
                                            
                                            auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1](InGap& obj) {return obj.getuId() == sId1;}); // given a uId, find it in gaps
                                            
                                            if (gId != inGaps->end()) {
                                            
                                                path.add('G', sId1, '0', start, end);
                                            
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
                                    
                                    inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                    
                                    arguments.push_back(newLine); // last column
                                    
                                    sId1Or = arguments[2][0]; // get sequence orientation in the edge
                                    
                                    seqHeader = arguments[1];
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash1(seqHeader, uId); // header to hash table
                                        inSequences.insertHash2(uId, seqHeader); // header to hash table
                                    
                                        sId1 = uId;
                                        
                                        uId++;
                                        
                                        inSequences.setuId(uId); // we have touched a segment need to increase the unique segment counter
                                        
                                    }else{
                                        
                                        sId1 = got->second;
                                        
                                    }
                                    
                                    sId2Or = arguments[4][0]; // get sequence orientation in the edge
                                    
                                    seqHeader = arguments[3];
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash.find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash.end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash1(seqHeader, uId); // header to hash table
                                        inSequences.insertHash2(uId, seqHeader); // header to hash table
                                    
                                        sId2 = uId;
                                        
                                        uId++;
                                        
                                        inSequences.setuId(uId); // we have touched a segment need to increase the unique segment counter
                                        
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
                                    inSequences.setuId(uId+1);
                                    
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
                                            
                                            inSequences.insertHash1(component, uId); // header to hash table
                                            inSequences.insertHash2(uId, component); // header to hash table
                                        
                                            sId1 = uId;
                                            
                                            inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                                            
                                        }else{
                                            
                                            sId1 = got->second;
                                            
                                        }
                                        
                                        std::vector<InSegment>* inSegments = inSequences.getInSegments();
                                        std::vector<InGap>* inGaps = inSequences.getInGaps();
                                        
                                        auto sId = find_if(inSegments->begin(), inSegments->end(), [sId1](InSegment& obj) {return obj.getuId() == sId1;}); // given a uId, find it in nodes
                                    
                                        if (sId != inSegments->end()) {
                                            
                                            path.add('S', sId1, sId1Or);
                                             
                                        }else{
                                            
                                            auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1](InGap& obj) {return obj.getuId() == sId1;}); // given a uId, find it in gaps
                                            
                                            if (gId != inGaps->end()) {
                                            
                                                path.add('G', sId1, '0');
                                            
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
        
        if (!instructions.empty()) {
            
            verbose("\nStarted instruction execution");
        
            SAK sak; // create a new swiss army knife
            
            for (Instruction instruction : instructions) { // execute swiss army knife instructions
                
                sak.executeInstruction(inSequences, instruction);
                
                verbose(instruction.action + " instruction executed");
                
            }
        
        }

        if (!iAgpFileArg.empty() || (isPipe && (pipeType == 'a'))) {
            
            std::string pHeaderNew, pHeader1, pHeader2, gHeader, instruction, coord1, coord2;
            char pId1Or, pId2Or;
            
            unsigned int pUId1 = 0, pUId2 = 0, gUId = 0, dist = 0, seqLen, pathLen, start1 = 0, end1 = 0, start2 = 0, end2 = 0;
            phmap::flat_hash_map<std::string, unsigned int> hash;
            phmap::flat_hash_map<std::string, unsigned int>::const_iterator got;
            
            std::vector<std::string> arguments; // line arguments
            std::vector<unsigned int> subsetted; // vector of paths flagged to be removed only at the end
            
            if (isPipe && (pipeType == 'a')) {
                
                stream = make_unique<std::istream>(std::cin.rdbuf());
                
            }else{
                
                stream = make_unique<std::ifstream>(std::ifstream(iAgpFileArg));
                
            }
            
            while (getline(*stream, line)) {
                
                std::istringstream iss(line); // line to string
                
                arguments = readDelimited(line, "\t", "#"); // read the columns in the line
                
                if (arguments.size() == 0) {continue;}
                
                pHeaderNew = arguments[0]; // this is the current header
                
                if (arguments[4] == "W") { // this is an old path
                    
                    pHeader1 = arguments[5];
                    pId1Or = arguments[8][0];
                    
                    hash = inSequences.getHash1();
                    
                    got = hash.find(pHeader1); // get the headers to uIds table
                    
                    if (got != hash.end()) { // this is not the first time we see this path
                        
                        pUId1 = got->second;
                        
                    }else{
                        
                        fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader1.c_str()); // sequence not found
                        
                        continue;
                        
                    }
                    
                    pathLen = inSequences.pathLength(pUId1);
                    
                    start1 = stoi(arguments[6]);
                    end1 = stoi(arguments[7]);
                    
                    seqLen = end1 - start1 + 1;
                    
                    if(seqLen != pathLen) {

                        fprintf(stderr, "Warning: sequence length (%u) differs from path length (%u). Subsetting (%s).\n", seqLen, pathLen, pHeader1.c_str());

                    }else{
                        
                        start1 = 0;
                        end1 = 0;
                        
                    }
                    
                    std::streampos oldpos = stream->tellg();  // stores the position before we read another line
                    
                    getline(*stream, line);
                    
                    arguments = readDelimited(line, "\t", "#"); // read the next sequence
                    
                    if(pHeaderNew != arguments[0]) { // if this path does not need to be joined to anything that follows, we simply rename it
                        
                        inSequences.renamePath(pUId1, pHeaderNew);
                        
                        if(pId1Or == '-') {
                            
                            inSequences.revComPath(pUId1);
                            
                        }
                        
                        if(seqLen != pathLen) { // if it also needs to be trimmed
                            
                            inSequences.trimPath(pUId1, start1, end1);
                            
                        }
                        
                    }else if (seqLen != pathLen){ // if the path needs to be joined to something and is subsetted
                        
//                        subsetted.push_back(pUId1);
                        
                    }
                    
                    stream->seekg(oldpos); // reset stream to previous line
                    
                }else if(arguments[4] == "N"){

                    hash = inSequences.getHash1();
                    
                    got = hash.find(pHeader1); // get the headers to uIds table (remove sequence orientation in the gap first)
                    
                    if (got == hash.end()) { // this is the first time we see this path
                        
                        fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader1.c_str()); // if the preceding sequence was not found we do not introduce a gap
                        
                        continue;
                        
                    }
                    
                    gUId = inSequences.getuId();
                    
                    if (arguments[6] == "scaffold") {
                    
                        gHeader = "gap"+std::to_string(gUId);
                    
                    }else{
                        
                        gHeader = arguments[6];
                        
                    }
                    
                    inSequences.insertHash1(gHeader, gUId); // header to hash table
                    inSequences.insertHash2(gUId, gHeader); // uId to hash table
                    
                    inSequences.setuId(gUId+1);
                    
                    dist = stoi(arguments[5]);
                    
                    getline(*stream, line);
                    
                    arguments = readDelimited(line, "\t", "#"); // read the next sequence
                    
                    if (arguments.size() == 0) {continue;}
                    
                    pHeader2 = arguments[5];
                    pId2Or = arguments[8][0];
                    
                    hash = inSequences.getHash1();
                    
                    got = hash.find(pHeader2); // get the headers to uIds table (remove sequence orientation in the gap first)
                    
                    if (got != hash.end()) { // this is not the first time we see this path
                        
                        pUId2 = got->second;
                        
                    }else{
                        
                        fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader2.c_str()); // sequence not found
                        
                        continue;
                        
                    }
                    
                    pathLen = inSequences.pathLength(pUId2);
                    
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
            
            if (subsetted.size() > 0) { // if sequences were subsetted, remove the original paths
                
                // avoid duplicates
                std::vector<unsigned int>::iterator ip;
                sort(subsetted.begin(), subsetted.end());
                ip = std::unique(subsetted.begin(), subsetted.end());
                subsetted.resize(std::distance(subsetted.begin(), ip));
                
                for (unsigned int pUId : subsetted) {
                    
                    inSequences.removePath(pUId);
                    
                }
                
            }
            
        }
            
        inSequences.updateScaffoldStats();
        
        verbose("Updated scaffold statistics");
            
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
            
        }

        return inSequences;
        
    }
        
    bool includeExcludeAppend(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality = NULL) {
 
        bedIncludeListHeaders = bedIncludeList.getSeqHeaders();
        bedExcludeListHeaders = bedExcludeList.getSeqHeaders();
        bool outSeq;
        
        verbose("Processing sequence: " + *seqHeader);
        
        if   (bedIncludeList.empty() &&
              bedExcludeList.empty()) {
            
            inSequences->appendSequence(seqHeader, seqComment, inSequence, inSequenceQuality);
            
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
            
            std::cout<<inSequence->size()<<std::endl;
                
            if (outSeq && inSequence->size()>0) {
                
                if (offset>0) {
                
                    inSequence->erase(offset, inSequence->size()-offset);
                    
                    if (inSequenceQuality != NULL) {
                    
                        inSequenceQuality->erase(offset, inSequenceQuality->size()-offset);
                        
                    }
                    
                }
                
                inSequences->appendSequence(seqHeader, seqComment, inSequence, inSequenceQuality);
            
            }else {
                
                verbose("Sequence entirely removed as a result of include: " + *seqHeader);
                
            }
                
        }else if(bedIncludeList.empty() &&
                !bedExcludeList.empty()) {
            
            pos = 0;
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
            
                inSequences->appendSequence(seqHeader, seqComment, inSequence, inSequenceQuality);
            
            }else {
                
                verbose("Sequence entirely removed as a result of exclude: " + *seqHeader);
                
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
                    
                    inSequences->appendSequence(seqHeader, seqComment, inSequence, inSequenceQuality);
                    
        }
        
        return false;
        
    }
    
    bool includeExcludeAppendSegment(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality = NULL) {
 
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
