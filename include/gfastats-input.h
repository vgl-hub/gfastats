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
                        
                        includeExcludeAppend(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList);
                        
                        getline(*stream, newLine);
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }
                        
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

                        includeExcludeAppend(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList, &inSequenceQuality);
                        
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
                    unsigned int sId1 = 0, sId2 = 0, dist = 0;
                    phmap::flat_hash_map<std::string, unsigned int> hash;
                    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got;
                    
                    unsigned int lineN = 1;
                    unsigned int uId = 0, guId = 0, euId = 0;
                    
                    std::string delimiter = "\t";
                    std::vector<std::string> arguments; // process the columns of each row
                    
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
                                    
                                    includeExcludeAppendSegment(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList);
                                    
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
                                    
                                    strtok(strdup(newLine.c_str()),"\t"); // process first line
                                    
                                    seqHeader = strtok(NULL,"\t");
                                    
                                    uId = inSequences.getuId();
                                    inSequences.setuId(uId+1);
                                    
                                    path.newPath(uId, seqHeader);
                                    
                                    s = strtok(NULL,"\t");
                                    
                                    delimiter = " ";
                                    
                                    arguments.clear();
                                    
                                    while ((pos = s.find(delimiter)) != std::string::npos) {
                                        
                                        arguments.push_back(s.substr(0, pos));
                                        
                                        s.erase(0, pos + delimiter.length());
                                    
                                    }
                                    
                                    arguments.push_back(s); // last column
                                    
                                    for (std::string component : arguments) {
                                        
                                        sId1Or = component.back(); // get sequence orientation
                                        
                                        component.pop_back();
                                    
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
                                            
                                            path.addToPath('S', sId1, sId1Or);
                                             
                                        }else{
                                            
                                            auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1](InGap& obj) {return obj.getuId() == sId1;}); // given a uId, find it in gaps
                                            
                                            if (gId != inGaps->end()) {
                                            
                                                path.addToPath('G', sId1);
                                            
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
                            
                            switch (newLine[0]) {
                                    
                                case 'S': {
                                    
                                    strtok(strdup(newLine.c_str()),"\t"); //process first line
                                    h = strtok(NULL,"\t");
                                    
                                    seqHeader = h;
                                    
                                    s = strtok(NULL,"\t");
                                    inSequence = s;
                                    
                                    c = strtok(NULL,"\t");
                                    if (c != NULL) {
                                        
                                        seqComment = std::string(c);
                                        
                                    }
                                    
                                    includeExcludeAppendSegment(&inSequences, &seqHeader, &seqComment, &inSequence, bedIncludeList, bedExcludeList);
                                    
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
                                        
                                        component.pop_back();
                                    
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
                                            
                                            path.addToPath('S', sId1, sId1Or);
                                             
                                        }else{
                                            
                                            auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1](InGap& obj) {return obj.getuId() == sId1;}); // given a uId, find it in gaps
                                            
                                            if (gId != inGaps->end()) {
                                            
                                                path.addToPath('G', sId1);
                                            
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
            
            if (instructions.empty() && inSequences.getGaps().size() > 0) { // it only makes sense to update the stats if we are not manipulating the sequence and there are gaps
                
                inSequences.updateScaffoldStats();
                
                verbose("Updated scaffold statistics");
                
            }
                
        }else{

            printf("Stream not successful: %s", iSeqFileArg.c_str());
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
            
            inSequences.updateScaffoldStats();
            
            verbose("Updated scaffold statistics");
        
        }

        if (!iAgpFileArg.empty() || (isPipe && (pipeType == 'a'))) {
            
            inSequences.clearPaths();
            verbose("Existing paths removed");
            
            inSequences.clearGaps();
            verbose("Existing gaps removed");
            
            std::string sHeader, prevsHeader = "", gHeader;
            char sId1Or, sId2Or;
            
            InGap gap;
            InPath path;
            unsigned int uId = 0, sId1 = 0, sId2 = 0, guId = 0, dist = 0;
            phmap::flat_hash_map<std::string, unsigned int> hash;
            phmap::flat_hash_map<std::string, unsigned int>::const_iterator got;
            
            std::vector<std::string> arguments; // line arguments
            
            if (isPipe && (pipeType == 'a')) {
                
                stream = make_unique<std::istream>(std::cin.rdbuf());
                
            }else{
                
                stream = make_unique<std::ifstream>(std::ifstream(iAgpFileArg));
                
            }
            
            while (getline(*stream, line)) {
                
                std::istringstream iss(line); // line to string
                
                arguments = readDelimited(line, "\t"); // read the columns in the line
                
                if (arguments[0] != prevsHeader) {
                    
                    if(prevsHeader != ""){
                    
                        inSequences.addPath(path);
                    
                    }
                        
                    uId = inSequences.getuId();
                    
                    inSequences.insertHash1(arguments[0], uId); // header to hash table
                    inSequences.insertHash2(uId, arguments[0]); // uid to hash table
                    
                    inSequences.setuId(uId+1);
                    
                    path.newPath(uId, arguments[0]); // set uid and header of the new path
                    
                    prevsHeader = arguments[0];
                    
                }

                if (arguments[4] == "W") {
                    
                    sHeader = arguments[5];
                    sId1Or = arguments[8][0];
                    
                    hash = inSequences.getHash1();
                    
                    got = hash.find(sHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                    
                    if (got != hash.end()) { // this is the first time we see this segment
                        
                        sId1 = got->second;
                        
                    }else{
                        
                        printf("Error: contig missing from the segment set (%s).\n", sHeader.c_str()); exit(1);
                        
                    }
                    
                    if(stoi(arguments[7]) != inSequences.getInSegment(sId1).getSegmentLength()) {
                        
                        printf("Error: contig length differs from segment length (%s).\n", sHeader.c_str()); exit(1);
                        
                    }
                    
                    path.addToPath('S', sId1, sId1Or);
                    
                }else if(arguments[4] == "N"){
                    
                    uId = inSequences.getuId();
                    
                    if (arguments[6] == "scaffold") {
                    
                        gHeader = "gap"+std::to_string(uId);
                    
                    }else{
                        
                        gHeader = arguments[6];
                        
                    }
                        
                    inSequences.insertHash1(gHeader, uId); // header to hash table
                    inSequences.insertHash2(uId, gHeader); // uId to hash table
                    
                    guId = uId; // since I am still reading segments I need to keep this fixed
                    
                    dist = stoi(arguments[5]);
                    
                    inSequences.setuId(uId+1); // we have touched a feature need to increase the unique feature counter
                    
                    std::streampos oldpos = stream->tellg();  // stores the position before we read another line
                    
                    getline(*stream, line);
                    
                    arguments = readDelimited(line, "\t"); // read the next sequence
                    
                    if (arguments[0] != prevsHeader) { // this path ends with a gap
                        
                        gap.newGap(guId, sId1, sId1, '+', '+', dist, gHeader); // now that we have all the information, create the gap
                        
                        inSequences.addGap(gap);
                        
                        path.addToPath('G', guId); // add the gap
                        
                        stream->seekg(oldpos); // reset stream to previous line
                    
                        continue;
                        
                    }
                    
                    sHeader = arguments[5];
                    sId2Or = arguments[8][0];
                    
                    hash = inSequences.getHash1();
                    
                    got = hash.find(sHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                    
                    if (got != hash.end()) { // this is the first time we see this segment
                        
                        sId2 = got->second;
                        
                    }else{
                        
                        printf("Error: contig missing from the segment set (%s).\n", sHeader.c_str()); exit(1);
                        
                    }
                    
                    if(stoi(arguments[7]) != inSequences.getInSegment(sId2).getSegmentLength()) {
                        
                        printf("Error: contig length differs from segment length (%s).\n", sHeader.c_str()); exit(1);
                        
                    }
                    
                    gap.newGap(guId, sId1, sId2, sId1Or, sId2Or, dist, gHeader); // now that we have all the information, create the gap
                    
                    inSequences.addGap(gap);
                    
                    path.addToPath('G', guId); // add the gap
                    
                    path.addToPath('S', sId2, sId2Or); // now add the segment it leads to
                    
                    sId1Or = sId2Or;
                    sId1 = sId2;
                    
                }
                
            }
            
            inSequences.addPath(path); // push the last path
            
            inSequences.updateScaffoldStats();
            
            verbose("Updated scaffold statistics");
            
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
            
        }

        return inSequences;
        
    }
        
    void includeExcludeAppend(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality = NULL) {
 
        bedIncludeListHeaders = bedIncludeList.getSeqHeaders();
        bedExcludeListHeaders = bedExcludeList.getSeqHeaders();
        bool outSeq;
        
        if   (bedIncludeList.empty() &&
              bedExcludeList.empty()) {
            
            verbose("Processing sequence: " + *seqHeader);
            
            inSequences->appendSequence(seqHeader, seqComment, inSequence, inSequenceQuality);
            
        }else if
            (!bedIncludeList.empty() &&
             bedExcludeList.empty()) {
            
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
                
                inSequences->appendSequence(seqHeader, seqComment, inSequence, inSequenceQuality);
            
            }else {
                
                verbose("Scaffold entirely removed as a result of include: " + *seqHeader);
                
            }
                
        }else if
            (bedIncludeList.empty() &&
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
                
                verbose("Scaffold entirely removed as a result of exclude: " + *seqHeader);
                
            }
                    
        }else if
                (!bedIncludeList.empty() &&
                 !bedExcludeList.empty() &&
                 std::find(bedIncludeListHeaders.begin(), bedIncludeListHeaders.end(), *seqHeader) != bedIncludeListHeaders.end() &&
                 std::find(bedExcludeListHeaders.begin(), bedExcludeListHeaders.end(), *seqHeader) == bedExcludeListHeaders.end()) {
                    
                    inSequences->appendSequence(seqHeader, seqComment, inSequence, inSequenceQuality);
                    
        }
        
    }
    
    void includeExcludeAppendSegment(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality = NULL) {
 
        bedIncludeListHeaders = bedIncludeList.getSeqHeaders();
        bedExcludeListHeaders = bedExcludeList.getSeqHeaders();
        bool outSeq;
        
        if   (bedIncludeList.empty() &&
              bedExcludeList.empty()) {
            
            inSequences->appendSegment(seqHeader, seqComment, inSequence, inSequenceQuality);
            
        }else if
            (!bedIncludeList.empty() &&
             bedExcludeList.empty()) {
            
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
                
        }else if
            (bedIncludeList.empty() &&
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
                    
                    inSequences->appendSegment(seqHeader, seqComment, inSequence, inSequenceQuality);
                    
        }
        
    }
    
};

#endif /* GFASTATS_INPUT_H */
