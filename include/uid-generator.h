#ifndef UID_GENERATOR_H
#define UID_GENERATOR_H

class UIdGenerator {
private:
    int prevUId;
    int currUId;
    int nextUId;

public:
    UIdGenerator();

    // returns the next uId without generating a new one
    int peek();

    // returns the current uId without generating a new uId
    int get();

    // returns the previously generated uId
    int prev();

    // returns the current uId and generates a new uId
    int next();
};

#endif /* UID_GENERATOR_H */
