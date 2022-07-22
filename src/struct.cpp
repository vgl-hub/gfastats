#include <stdlib.h>
#include <string>
#include <vector>

#include "bed.h"
#include "struct.h"

std::string UserInput::file(char type) {
    
    std::string filename;
    
    switch (type) {
        case 'f':
            filename = iSeqFileArg;
            break;
        case 'r':
            filename = iReadFileArg;
            break;
        case 'i':
            filename = iBedIncludeFileArg;
            break;
        case 'e':
            filename = iBedExcludeFileArg;
            break;
        case 'a':
            filename = iAgpFileArg;
            break;
        case 'k':
            filename = iSakFileArg;
            break;
    }
    
    return filename;
    
}

Sequences::~Sequences()
{
    for (Sequence* p : sequences)
        delete p;
    
}

Sequence::~Sequence()
{
    delete sequence;
    delete sequenceQuality;
}

bool Edge::operator==(const Edge& e) const {
    return orientation0 == e.orientation0 && orientation1 == e.orientation1 && id == e.id;
}
