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
    
    void readFiles(InSequences &inSequences, std::string &iSeqFileArg, std::string &iSakFileArg, std::string &iAgpFileArg, std::string &iBedIncludeFileArg, std::string &iBedExcludeFileArg, BedCoordinates &bedIncludeList, bool isPipe, char &pipeType, std::string sortType, std::string &iReadFileArg);
    
    Sequence includeExcludeSeq(std::string seqHeader, std::string seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality = NULL);

    Sequence includeExcludeSeg(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality = NULL, std::vector<Tag>* inTags = NULL);
    
};

#endif /* GFASTATS_INPUT_H */
