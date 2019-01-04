//
// Created by coki on 04.01.19..
//

#include "Sequence.h"

const std::string &Sequence::getDescription() const {
    return description;
}

const std::vector<char> &Sequence::getSequence() const {
    return sequence;
}

void Sequence::setSequence(const std::vector<char>& s) {
    sequence = s;
}

void Sequence::setDescription(const std::string &d) {
    description = d;
}
