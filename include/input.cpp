//
//  gfastats-input.h
//  gfastats
//
//  Created by Giulio Formenti on 1/16/22.
//

#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>

#include <iostream>
#include <fstream>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "gfastats-struct.h"
#include "gfastats-functions.h" // global functions

#include "gfastats-log.h"
#include "gfastats-global.h"
#include "uid-generator.h"

#include "gfa-lines.h"
#include "reads.h"

#include "threadpool.h"
#include "gfa.h"
#include "gfastats-sak.h" // swiss army knife

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "stream-obj.h"

#include "input-agp.h"
#include "input-filters.h"
#include "input.h"

void Input::load(UserInput userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(InReads& inReads) {
    
    if (!userInput.iReadFileArg.empty()) {

        inReads.load(userInput);

    }

}
    
void Input::read(InSequences& inSequences) {

    if (!userInput.iSakFileArg.empty() || userInput.pipeType == 'k') {
        
        stream = streamObj.openStream(userInput, 'k');
        
        SAK sak; // create a new swiss army knife
        
        while (getline(*stream, line)) {
            
            std::istringstream iss(line);
            
            instructions.push_back(sak.readInstruction(line)); // use the swiss army knife to read the instruction
            
        }
        
        lg.verbose("Finished reading SAK instructions");
        
    }
    
    if (!userInput.iBedIncludeFileArg.empty() || userInput.pipeType == 'i') {
        
        stream = streamObj.openStream(userInput, 'i');
        
        while (getline(*stream, line)) {
            
            std::istringstream iss(line);
            iss >> bedHeader >> begin >> end;
            
            userInput.bedIncludeList.pushCoordinates(bedHeader, begin, end);
            begin = 0, end = 0;
            
        }
        
        lg.verbose("Finished reading BED include list");
        
    }
    
    BedCoordinates bedExcludeList;
    
    if (!userInput.iBedExcludeFileArg.empty() || userInput.pipeType == 'e') {
        
        stream = streamObj.openStream(userInput, 'e');
        
        while (getline(*stream, line)) {
            
            std::istringstream iss(line);
            iss >> bedHeader >> begin >> end;
            
            bedExcludeList.pushCoordinates(bedHeader, begin, end);
            begin = 0, end = 0;
            
        }
        
        lg.verbose("Finished reading BED exclude list");
        
    }
    
    if (!userInput.iSeqFileArg.empty() || userInput.pipeType == 'f') {
        
        stream = streamObj.openStream(userInput, 'f');
        
        lg.verbose("Created stream object for input assembly file");
        lg.verbose("Detected stream type (" + streamObj.type() + ").\nStreaming started.");
        
        if (stream) {
            
            switch (stream->peek()) {
                    
                case '>': {
                    
                    stream->get();
                    
                    while (getline(*stream, newLine)) {
                        
                        if(userInput.bedIncludeList.size() - bedExcludeList.size() != 0 && userInput.bedIncludeList.size() - bedExcludeList.size() == inSequences.getPathN()) { // we have all the sequences needed
                            lg.verbose("Found all sequences, stop streaming input");
                            break;
                        
                        }
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }
                        
                        std::string* inSequence = new std::string;
                        
                        getline(*stream, *inSequence, '>');
                        
                        lg.verbose("Individual fasta sequence read");
                        
                        Sequence* sequence = includeExcludeSeq(seqHeader, seqComment, inSequence, userInput.bedIncludeList, bedExcludeList);
                        
                        if (sequence != NULL) {
                            
                            sequence->seqPos = seqPos; // remember the order
                            
                            inSequences.appendSequence(sequence);
                            
                            seqPos++;
                            
                        }
                        
                    }
                    
                    break;
                }
                case '@': {
                    
                    while (getline(*stream, newLine)) { // file input
                        
                        if(userInput.bedIncludeList.size() - bedExcludeList.size() != 0 && userInput.bedIncludeList.size() - bedExcludeList.size() == inSequences.getPathN()) { // we have all the sequences needed
                            lg.verbose("Found all sequences, stop streaming input");
                            break;
                        
                        }
                        
                        newLine.erase(0, 1);
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }
                        
                        std::string* inSequence = new std::string;
                        getline(*stream, *inSequence);
                        
                        getline(*stream, newLine);
                        
                        std::string* inSequenceQuality = new std::string;
                        getline(*stream, *inSequenceQuality);

                        Sequence* sequence = includeExcludeSeq(seqHeader, seqComment, inSequence, userInput.bedIncludeList, bedExcludeList, inSequenceQuality);
                        
                        if (sequence != NULL) {
                            
                            sequence->seqPos = seqPos; // remember the order
                        
                            inSequences.appendSequence(sequence);
                            
                            seqPos++;
                            
                        }
                        
                    }
                    
                    break;
                    
                }
                default: {
                    
                    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
                    
                    std::string h_col1, h_col2, h_col3, s, version, gHeader, eHeader, cigar, startS, endS;
                    char sId1Or, sId2Or;
                    
                    InGap gap;
                    InEdge edge;
                    InPath path;
                    unsigned int sId1 = 0, sId2 = 0, dist = 0, start = 0, end = 0, gapN = 0;
                    phmap::flat_hash_map<std::string, unsigned int>* hash;
                    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got;
                    
                    unsigned int uId = 0, guId = 0, euId = 0;
                    
                    bool isSegment = false;
                    
                    std::vector<std::string> arguments, components, tagValues; // process the columns of each row
                    
                    std::vector<Tag> inTags;
                    Tag tag;
                    
                    getline(*stream, newLine);
                    
                    firstLine = newLine;
                    
                    arguments = readDelimited(newLine, "\t");
                    
                    if (arguments[0] == "H") {
                        
                        h_col2 = arguments[1]; // read header col2
                        
                        arguments = readDelimited(newLine, ":");
                        
                        if (arguments[2] != "") {
                            
                            version = arguments[2];
                            lg.verbose("GFA version: " + version);
                            
                        }else{
                            
                            lg.verbose("Failed to parse GFA version. Please check your header.");
                            
                        }
                    
                    }else{
                            
                        lg.verbose("Cannot recognize GFA version from first line. Trying to detect from content.");
                        
                        if (arguments[0] == "S") {
                            
                            if (isInt(arguments[2]) || (arguments[2] == "*" && arguments[3].find(":") == std::string::npos)) {
                                
                                version = '2';
                                lg.verbose("Proposed GFA version: " + version);
                                
                            }else{
                                
                                version = '1';
                                lg.verbose("Proposed GFA version: " + version);
                                
                            }
                            
                        }else if (arguments[0] == "G" || arguments[0] == "O") {
                            
                            version = '2';
                            lg.verbose("Proposed GFA version: " + version);
                            
                        }
                        
                        stream->seekg(0); // return to begin of file after detection
                            
                    }
                    
                    if (version[0] == '2') { // GFA2
                    
                        while (getline(*stream, newLine)) {
                            
                            if (stopStream) {break;}
                            
                            switch (newLine[0]) {
                                    
                                case 'S': {
                                    
                                    arguments = readDelimited(newLine, "\t");
                                    
                                    seqHeader = arguments[1];
                                    
                                    std::string* inSequence = new std::string;
                                    
                                    *inSequence = arguments[3];
                                    
                                    inTags.clear();
                                    
                                    for (unsigned int i = 4; i < arguments.size(); i++) {
                                        
                                        tagValues = readDelimited(arguments[i], ":");
                                        
                                        tag.label[0] = tagValues[0][0];
                                        tag.label[1] = tagValues[0][1];
                                        tag.type = tagValues[1][0];
                                        tag.content = tagValues[2];
                                    
                                        inTags.push_back(tag);
                                    
                                    }
                                    
                                    Sequence* sequence = includeExcludeSeg(&inSequences, &seqHeader, &seqComment, inSequence, userInput.bedIncludeList, bedExcludeList);
                                    
                                    if (sequence != NULL) {
                                        
                                        sequence->seqPos = seqPos; // remember the order
                                    
                                        inSequences.appendSegment(sequence, inTags);
                                        seqPos++;
                                        
                                    }
                                    
                                    break;
                                    
                                }
                                case 'G': {
                                    
                                    lck.lock();
                                    
                                    if(verbose_flag) {std::cerr<<"\n\n";};
                                    
                                    arguments = readDelimited(newLine, "\t");
                                    
                                    gHeader = arguments[1];
                                    
                                    uId = inSequences.getuId();
                                    
                                    inSequences.insertHash(gHeader, uId);
                                    
                                    guId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                    
                                    sId1Or = arguments[2].back(); // get sequence orientation in the gap
                                    
                                    seqHeader = std::string(arguments[2]);
                                    seqHeader.pop_back();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash->end()) { // this is the first time we see this segment
                                        
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
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash->end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId2 = uId;
                                        
                                        inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                        
                                    }else{
                                        
                                        sId2 = got->second;
                                        
                                    }
                                    
                                    dist = stoi(arguments[4]);
                                    
                                    lg.verbose("Processing gap " + gHeader + " (uId: " + std::to_string(uId) + ")");
                                    
                                    inTags.clear();
                                    
                                    for (unsigned int i = 5; i < arguments.size(); i++) {
                                        
                                        tagValues = readDelimited(arguments[i], ":");
                                        
                                        tag.label[0] = tagValues[0][0];
                                        tag.label[1] = tagValues[0][1];
                                        tag.type = tagValues[1][0];
                                        tag.content = tagValues[2];
                                    
                                        inTags.push_back(tag);
                                    
                                    }
                                    
                                    gap.newGap(guId, sId1, sId2, sId1Or, sId2Or, dist, gHeader, inTags);
                                    
                                    inSequences.addGap(gap);
                                    
                                    lck.unlock();
                                                 
                                    break;
                                    
                                }

                                case 'E': {
                                    
                                    lck.lock();
                                    
                                    if(verbose_flag) {std::cerr<<"\n\n";};
                                    
                                    arguments = readDelimited(newLine, "\t");
                                    
                                    eHeader = arguments[1];
                                    
                                    uId = inSequences.getuId();
                                    
                                    inSequences.insertHash(eHeader, uId);
                                    
                                    euId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                    
                                    sId1Or = arguments[2].back(); // get sequence orientation in the edge
                                    
                                    seqHeader = std::string(arguments[2]);
                                    seqHeader.pop_back();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash->end()) { // this is the first time we see this segment
                                        
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
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash->end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId2 = uId;
                                        
                                        inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                        
                                    }else{
                                        
                                        sId2 = got->second;
                                        
                                    }
      
                                    cigar = arguments[8];
                                    
                                    inTags.clear();
                                    
                                    for (unsigned int i = 9; i < arguments.size(); i++) {
                                        
                                        tagValues = readDelimited(arguments[i], ":");
                                        
                                        tag.label[0] = tagValues[0][0];
                                        tag.label[1] = tagValues[0][1];
                                        tag.type = tagValues[1][0];
                                        tag.content = tagValues[2];
                                    
                                        inTags.push_back(tag);
                                    
                                    }
                                    
                                    edge.newEdge(euId, sId1, sId2, sId1Or, sId2Or, cigar, eHeader, inTags);
                                    
                                    inSequences.appendEdge(edge);
                                    
                                    lck.unlock();
                                                 
                                    break;
                                    
                                }
                                    
                                case 'O': {
                                    
                                    lck.lock();
                                    
                                    if(verbose_flag) {std::cerr<<"\n\n";};
                                    
                                    arguments = readDelimited(newLine, "\t");
                                    
                                    seqHeader = arguments[1];
                                    
                                    uId = inSequences.getuId();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table to look for the header
                                    
                                    if (got == hash->end()) { // this is the first time we see this header
                                        
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
                                            isSegment = true;
                                            
                                        }else{
                                            
                                            isSegment = false;
                                            
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
                                        
                                            lg.verbose("Adding only coordinates " + std::to_string(start) + ":" + std::to_string(end) + "(" + component + ")");
                                            
                                        }
                                    
                                        hash = inSequences.getHash1();
                                        
                                        got = hash->find(component); // get the headers to uIds table (remove sequence orientation in the gap first)
                                        
                                        if (got == hash->end()) { // this is the first time we see this component
                                            
                                            uId = inSequences.getuId();
                                            
                                            inSequences.insertHash(component, uId);
                                        
                                            sId1 = uId;
                                            
                                            inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                            
                                        }else{
                                            
                                            sId1 = got->second;
                                            
                                        }
                                    
                                        if (isSegment) {
                                            
                                            path.add(SEGMENT, sId1, sId1Or, start, end);
                                             
                                        }else{
                                            
                                            path.add(GAP, sId1, '0', start, end);
                                            
                                        }
                                        
                                    }
                                    
                                    for (unsigned int i = 2; i < arguments.size(); i++) {
                                        
                                        if (arguments[i].substr(0,3) == "C:Z") {
                                            
                                            seqComment = arguments[i];
                                            
                                            seqComment.erase(0,4);
                                            
                                            path.setComment(seqComment);
                                            
                                        }
                                        
                                    }
                                    
                                    inSequences.addPath(path);
                                    
                                    lck.unlock();
                                    
                                    break;
                                    
                                }
                                    
                            }
                            
                        }
                        
                    }else if (version[0] == '1') {
                    
                        while (getline(*stream, newLine)) {
                            
                            if (stopStream) {break;}
                            
                            switch (newLine[0]) {
                                    
                                case 'S': {
                                    
                                    arguments = readDelimited(newLine, "\t");
                                    
                                    seqHeader = arguments[1];
                                    
                                    std::string* inSequence = new std::string;
                                    
                                    *inSequence = arguments[2];
                                    
                                    inTags.clear();
                                    
                                    for (unsigned int i = 3; i < arguments.size(); i++) {
                                        
                                        tagValues = readDelimited(arguments[i], ":");
                                        
                                        tag.label[0] = tagValues[0][0];
                                        tag.label[1] = tagValues[0][1];
                                        tag.type = tagValues[1][0];
                                        tag.content = tagValues[2];
                                    
                                        inTags.push_back(tag);
                                    
                                    }
                                    
                                    Sequence* sequence = includeExcludeSeg(&inSequences, &seqHeader, &seqComment, inSequence, userInput.bedIncludeList, bedExcludeList, NULL, &inTags);
                                    
                                    if (sequence != NULL) {
                                        
                                        sequence->seqPos = seqPos; // remember the order
                                    
                                        inSequences.appendSegment(sequence, inTags);
                                        seqPos++;
                                        
                                    }
                                    
                                    break;
                                    
                                }
                                    
                                case 'J': {
                                    
                                    lck.lock();
                                    
                                    if(verbose_flag) {std::cerr<<"\n\n";};
                                    
                                    arguments = readDelimited(newLine, "\t");
                                    
                                    gHeader = "gap" + std::to_string(gapN);
                                    
                                    gapN++;
                                    
                                    uId = inSequences.getuId();
                                    
                                    inSequences.insertHash(gHeader, uId);
                                    
                                    guId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                    
                                    seqHeader = std::string(arguments[1]); // first component
                                    
                                    sId1Or = arguments[2][0]; // get orientation in the gap
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash->end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId1 = uId;
                                        
                                        inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                        
                                    }else{
                                        
                                        sId1 = got->second;
                                        
                                    }
                                    
                                    seqHeader = arguments[3]; // second component
                                    
                                    sId2Or = arguments[4][0]; // get orientation in the gap
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table (remove sequence orientation in the gap first)
                                    
                                    if (got == hash->end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId2 = uId;
                                        
                                        inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                        
                                    }else{
                                        
                                        sId2 = got->second;
                                        
                                    }
                                    
                                    dist = stoi(arguments[5]);
                                    
                                    lg.verbose("Processing gap " + gHeader + " (uId: " + std::to_string(uId) + ")");
                                    
                                    inTags.clear();
                                    
                                    for (unsigned int i = 6; i < arguments.size(); i++) {
                                        
                                        tagValues = readDelimited(arguments[i], ":");
                                        
                                        tag.label[0] = tagValues[0][0];
                                        tag.label[1] = tagValues[0][1];
                                        tag.type = tagValues[1][0];
                                        tag.content = tagValues[2];
                                    
                                        inTags.push_back(tag);
                                    
                                    }
                                    
                                    gap.newGap(guId, sId1, sId2, sId1Or, sId2Or, dist, gHeader, inTags);
                                    
                                    inSequences.addGap(gap);
                                    
                                    lck.unlock();
                                                 
                                    break;
                                    
                                }

                                case 'L': {
                                    
                                    lck.lock();
                                    
                                    if(verbose_flag) {std::cerr<<"\n\n";};
                                    
                                    arguments = readDelimited(newLine, "\t");

                                    uId = inSequences.getuId();
                                    
                                    euId = uId; // since I am still reading segments I need to keep this fixed
                                    
                                    inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                    
                                    sId1Or = arguments[2][0]; // get sequence orientation in the edge
                                    
                                    seqHeader = arguments[1];
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash->end()) { // this is the first time we see this segment
                                        
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
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table (remove sequence orientation in the edge first)
                                    
                                    if (got == hash->end()) { // this is the first time we see this segment
                                        
                                        uId = inSequences.getuId();
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                    
                                        sId2 = uId;
                                        
                                        uId++;
                                        
                                        inSequences.uId.next(); // we have touched a segment need to increase the unique segment counter
                                        
                                    }else{
                                        
                                        sId2 = got->second;
                                        
                                    }
                                    
                                    cigar = arguments[5];
                                    
                                    inTags.clear();
                                    
                                    for (unsigned int i = 6; i < arguments.size(); i++) {
                                        
                                        tagValues = readDelimited(arguments[i], ":");
                                        
                                        tag.label[0] = tagValues[0][0];
                                        tag.label[1] = tagValues[0][1];
                                        tag.type = tagValues[1][0];
                                        tag.content = tagValues[2];
                                    
                                        inTags.push_back(tag);
                                    
                                    }
                                    
                                    edge.newEdge(euId, sId1, sId2, sId1Or, sId2Or, cigar, "", inTags);
                                    
                                    inSequences.appendEdge(edge);
                                    
                                    lck.unlock();
                                                 
                                    break;
                                    
                                }

                                case 'P': {
                                    
                                    lck.lock();
                                    
                                    if(verbose_flag) {std::cerr<<"\n\n";};
                                    
                                    arguments = readDelimited(newLine, "\t");
                                    
                                    if (!(arguments[2].find(",") == std::string::npos)) break; // we are not reading edge paths yet
                                    
                                    seqHeader = arguments[1];
                                    
                                    uId = inSequences.getuId();
                                    
                                    hash = inSequences.getHash1();
                                    
                                    got = hash->find(seqHeader); // get the headers to uIds table to look for the header
                                    
                                    if (got == hash->end()) { // this is the first time we see this header
                                        
                                        inSequences.insertHash(seqHeader, uId);
                                        
                                    }else{
                                        
                                        fprintf(stderr, "Error: path name already exists (%s). Terminating.\n", seqHeader.c_str()); exit(1);
                                        
                                    }
                                    
                                    path.newPath(uId, seqHeader);
                                    
                                    inSequences.uId.next();
                                    
                                    components = readDelimited(arguments[2], ";");
                                    
                                    for (auto it = std::begin(components); it != std::end(components); ++it) {
                                        
                                        std::string component;
                                        
                                        if(it == std::begin(components) && *it == "") { // handle starting/ending gap
                                                
                                                component = *(std::next(it, 1));
                                                
                                                component.pop_back();
                                                
                                                got = hash->find(component); // get the headers to uIds table (remove sequence orientation in the gap first)
                                                
                                                if (got == hash->end()) { // this is the first time we see this segment
                                                    
                                                    fprintf(stderr, "Error1: cannot find next component in path (%s). Terminating.\n", component.c_str()); exit(1);
                                                    
                                                }else{
                                                    
                                                    sId1 = got->second;
                                                    
                                                }
                                                
                                                std::vector<InGap>* inGaps = inSequences.getInGaps();
                                                
                                                auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1](InGap& obj) {return obj.getsId1() == sId1 && obj.getsId2() == sId1;}); // given a uId, find it in gaps
                                                
                                                if (gId != inGaps->end()) {
                                                    
                                                    path.add(GAP, gId->getuId(), '0', start, end);
                                                    
                                                    lg.verbose("Adding gap to path with id:" + std::to_string(gId->getuId()));
                                                
                                                }
                                            
                                            ++it;
                                            
                                        }
                                        
                                        component = *it;
                                        
                                        sId1Or = component.back(); // get sequence orientation
                                        
                                        component.pop_back();
                                        
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
                                        
                                            lg.verbose("Adding only coordinates " + std::to_string(start) + ":" + std::to_string(end) + "(" + component + ")");
                                            
                                        }
                                    
                                        hash = inSequences.getHash1();
                                        
                                        got = hash->find(component); // get the headers to uIds table (remove sequence orientation in the gap first)
                                        
                                        if (got == hash->end()) { // this is the first time we see this segment
                                            
                                            uId = inSequences.getuId();
                                            
                                            inSequences.insertHash(component, uId);
                                        
                                            sId1 = uId;
                                            
                                            inSequences.uId.next(); // we have touched a feature need to increase the unique feature counter
                                            
                                        }else{
                                            
                                            sId1 = got->second;
                                            
                                        }
                                        
                                        std::vector<InGap>* inGaps = inSequences.getInGaps();
                                            
                                        path.add(SEGMENT, sId1, sId1Or, start, end);
                                        
                                        if (std::next(it, 1) != std::end(components)){
                                            
                                            component = *(std::next(it, 1));
                                            
                                            if (component != "") {
                                            
                                                component.pop_back();
                                                
                                                got = hash->find(component); // get the headers to uIds table (remove sequence orientation in the gap first)
                                                
                                                if (got == hash->end()) { // this is the first time we see this segment
                                                    
                                                    fprintf(stderr, "Error: cannot find next component in path (%s). Terminating.\n", component.c_str()); exit(1);
                                                    
                                                }else{
                                                    
                                                    sId2 = got->second;
                                                    
                                                }
                                                
                                                auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1,sId2](InGap& obj) {return (obj.getsId1() == sId1 && obj.getsId2() == sId2) || (obj.getsId1() == sId2 && obj.getsId2() == sId1);}); // given a uId, find it in gaps
                                                
                                                if (gId != inGaps->end()) {
                                                    
                                                    path.add(GAP, gId->getuId(), '0', start, end);
                                                    
                                                    lg.verbose("Adding gap to path with id:" + std::to_string(gId->getuId()));
                                                
                                                }
                                                
                                            }else{
                                                
                                                auto gId = find_if(inGaps->begin(), inGaps->end(), [sId1](InGap& obj) {return obj.getsId1() == sId1 && obj.getsId2() == sId1;}); // terminal gap
                                                
                                                if (gId != inGaps->end()) {
                                                    
                                                    path.add(GAP, gId->getuId(), '0', start, end);
                                                    
                                                    lg.verbose("Adding gap to path with id:" + std::to_string(gId->getuId()));
                                                
                                                }
                                                
                                            }
                                            
                                        }
                                        
                                    }
                                    
                                    for (unsigned int i = 2; i < arguments.size(); i++) {
                                        
                                        if (arguments[i].substr(0,3) == "C:Z") {
                                            
                                            seqComment = arguments[i];
                                            
                                            seqComment.erase(0,4);
                                            
                                            path.setComment(seqComment);
                                            
                                        }
                                        
                                    }
                                    
                                    inSequences.addPath(path);
                                    
                                    lck.unlock();
                                    
                                    break;
                                    
                                }
                                    
                            }
                            
                        }
                        
                        break;
                        
                    }
                    
                }
                
            }
            
            lg.verbose("End of file");
                
        }else{

            fprintf(stderr, "Stream not successful: %s", userInput.iSeqFileArg.c_str());
            exit(1);

        }
        
    }
    
    while (true) {
        
        if (threadPool.empty()) {threadPool.join(); break;}
        lg.verbose("Remaining jobs: " + std::to_string(threadPool.queueSize()), true);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        
    }
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    std::vector<Log> logs = inSequences.getLogs();
    
    //consolidate log
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    inSequences.sortSegmentsByOriginal();
    
    if (rmGaps_flag) {
     
        inSequences.removeTerminalGaps();
        
    }
    
    if (discoverPaths_flag) {
        
        inSequences.discoverPaths();
        
    }
    
    if (!instructions.empty()) {
        
        lg.verbose("\nStarted instruction execution");
    
        SAK sak; // create a new swiss army knife
        
        for (Instruction instruction : instructions) { // execute swiss army knife instructions
            
            sak.executeInstruction(inSequences, instruction);
            
            lg.verbose(instruction.action + " instruction executed");
            
        }
    
    }
    
    if (userInput.sortType == "ascending") {
        
        inSequences.sortPathsByNameAscending();
        
    }else if (userInput.sortType == "descending") {
        
        inSequences.sortPathsByNameDescending();
        
    }else if (userInput.sortType == "largest") {
        
        inSequences.sortPathsBySize(0);

    }else if (userInput.sortType == "smallest") {
        
        inSequences.sortPathsBySize(1);
        
    }else if (userInput.sortType != "none" && ifFileExists(userInput.sortType.c_str())){
            
        stream = std::make_unique<std::ifstream>(std::ifstream(userInput.sortType));
        
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

    if (!userInput.iAgpFileArg.empty() || userInput.pipeType == 'a') {
        
        readAgp(inSequences, userInput);
        
    }
        
    inSequences.updateStats();
    
}
