#ifndef INPUT_GFA_H
#define INPUT_GFA_H

void readGFA(InSequences& inSequences, UserInput& userInput, std::shared_ptr<std::istream> stream, BedCoordinates* bedExcludeList = NULL);

#endif /* INPUT_GFA_H */
