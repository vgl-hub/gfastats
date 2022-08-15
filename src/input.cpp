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
#include "struct.h"
#include "functions.h"

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "gfa-lines.h"

#include "threadpool.h"
#include "gfa.h"
#include "sak.h"

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "stream-obj.h"

#include "input-agp.h"
#include "input-filters.h"
#include "input-gfa.h"
#include "input.h"

void Input::load(UserInput userInput) {
    
    this->userInput = userInput;
    
}
    
void Input::read(InSequences& inSequences) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    threadPool.init(maxThreads); // initialize threadpool

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
                    
                    readGFA(inSequences, userInput, stream, &bedExcludeList);
                    
                }
                
            }
            
            lg.verbose("End of file");
                
        }else{

            fprintf(stderr, "Stream not successful: %s", userInput.iSeqFileArg.c_str());
            exit(1);

        }
        
    }
    
    jobWait(threadPool);
    threadPool.join();
    
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
    
    if (extractContigs_flag) {
        
        inSequences.clearGaps();
        inSequences.clearPaths();
        
    }
    
    if (extractContigs_flag || discoverPaths_flag) {
        
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
