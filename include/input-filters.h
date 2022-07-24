#ifndef INPUT_FILTERS_H
#define INPUT_FILTERS_H

Sequence* includeExcludeSeq(std::string seqHeader, std::string seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality = NULL);

Sequence* includeExcludeSeg(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates* bedExcludeList = NULL, std::string* inSequenceQuality = NULL);

#endif /* INPUT_FILTERS_H */
