#include "uid-generator.h"

UIdGenerator::UIdGenerator() {
    prevUId = 0;
    currUId = 0;
    nextUId = 1;
}

// returns the next uId without generating a new one
int UIdGenerator::peek() {
    return nextUId;
}

// returns the current uId without generating a new uId
int UIdGenerator::get() {
    return currUId;
}

// returns the previously generated uId
int UIdGenerator::prev() {
    return prevUId;
}

// returns the current uId and generates a new uId
int UIdGenerator::next() {
    prevUId = currUId;
    currUId = nextUId;
    nextUId ++;
    return prevUId;
}
