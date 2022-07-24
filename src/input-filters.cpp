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
#include "functions.h" // global functions

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "gfa-lines.h"

#include "threadpool.h"
#include "gfa.h"
#include "sak.h" // swiss army knife

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "stream-obj.h"

#include "input-agp.h"
#include "input-filters.h"

Sequence* includeExcludeSeq(std::string seqHeader, std::string seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality) {
    
    std::vector<std::string> bedIncludeListHeaders;
    std::vector<std::string> bedExcludeListHeaders;
    unsigned int pos = 0, cBegin = 0, cEnd = 0, offset = 0, prevCEnd = 0;

    bedIncludeListHeaders = bedIncludeList.getSeqHeaders();
    bedExcludeListHeaders = bedExcludeList.getSeqHeaders();
    bool outSeq = false;
    
    lg.verbose("Processing sequence: " + seqHeader);
    
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
                
    }else if
            (!bedIncludeList.empty() &&
             !bedExcludeList.empty() &&
             std::find(bedIncludeListHeaders.begin(), bedIncludeListHeaders.end(), seqHeader) != bedIncludeListHeaders.end() &&
             std::find(bedExcludeListHeaders.begin(), bedExcludeListHeaders.end(), seqHeader) == bedExcludeListHeaders.end()) {
                
                outSeq = true;
                
    }
    
    if (outSeq && inSequence->size()>0) {
    
        return new Sequence {seqHeader, seqComment, inSequence, inSequenceQuality};
    
    }else {
        
        lg.verbose("Sequence entirely removed as a result of BED filter: " + seqHeader);
        
        return NULL;
        
    }
    
}

Sequence* includeExcludeSeg(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates* bedExcludeList, std::string* inSequenceQuality) {
    
    std::vector<std::string> bedIncludeListHeaders;
    std::vector<std::string> bedExcludeListHeaders;
    unsigned int pos = 0, cBegin = 0, cEnd = 0, offset = 0, prevCEnd = 0;

    bedIncludeListHeaders = bedIncludeList.getSeqHeaders();
    
    if(bedExcludeList != NULL) {
    
        bedExcludeListHeaders = bedExcludeList->getSeqHeaders();
    
    }
        
    bool outSeq = false;
    
    lg.verbose("Processing sequence: " + *seqHeader);
    
    if   (bedIncludeList.empty() &&
          bedExcludeList->empty()) {
        
        outSeq = true;
        
    }else if(!bedIncludeList.empty() &&
              bedExcludeList != NULL &&
              bedExcludeList->empty()) {
        
        if(inSequences->getInSegments()->size() == bedIncludeList.size()) { // check if we retrieved all we needed
            
            lg.verbose("Found all sequences, stop streaming input");
            
            outSeq = true;
            
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
            
            outSeq = true;
        
        }else {
            
            lg.verbose("Sequence entirely removed as a result of include: " + *seqHeader);
            
        }
            
    }else if(bedIncludeList.empty() &&
             bedExcludeList != NULL &&
            !bedExcludeList->empty()) {
            
        offset = 0;
        outSeq = true;
        
        auto it = begin(bedExcludeListHeaders);
        
        while (it != end(bedExcludeListHeaders)) {
            
            it = std::find(it, bedExcludeListHeaders.end(), *seqHeader);
            
            if (it == end(bedExcludeListHeaders)) {
                
                break;
                
            }

            cBegin = bedExcludeList->getcBegin(pos);
            cEnd = bedExcludeList->getcEnd(pos);
            
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
        
            outSeq = true;
        
        }else {
            
            lg.verbose("Sequence entirely removed as a result of exclude: " + *seqHeader);
            
        }
                
    }else if
            (!bedIncludeList.empty() &&
              bedExcludeList != NULL &&
             !bedExcludeList->empty() &&
             std::find(bedIncludeListHeaders.begin(), bedIncludeListHeaders.end(), *seqHeader) != bedIncludeListHeaders.end() &&
             std::find(bedExcludeListHeaders.begin(), bedExcludeListHeaders.end(), *seqHeader) == bedExcludeListHeaders.end()) {
                
                if(inSequences->getInSegments()->size() == bedIncludeList.size()) { // check if we retrieved all we needed
                    
                    lg.verbose("Found all sequences, stop streaming input");
                    
                }
                
    }
    
    if (outSeq && inSequence->size()>0) {
    
        return new Sequence {*seqHeader, seqComment != NULL ? *seqComment : "", inSequence};
    
    }else {
        
        lg.verbose("Sequence entirely removed as a result of BED filter: " + *seqHeader);
        
        return NULL;
        
    }
    
}
