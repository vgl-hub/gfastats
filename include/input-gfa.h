#ifndef INPUT_GFA_H
#define INPUT_GFA_H

void readGFA(InSequences& inSequences, UserInput& userInput, BedCoordinates bedExcludeList, std::shared_ptr<std::istream> stream);

#endif /* INPUT_GFA_H */
