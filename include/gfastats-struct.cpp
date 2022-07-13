#include <string>
#include "gfastats-struct.h"

Sequence::~Sequence()
{
    delete sequence;
    delete sequenceQuality;
}

bool Edge::operator==(const Edge& e) const {
    return orientation0 == e.orientation0 && orientation1 == e.orientation1 && id == e.id;
}
