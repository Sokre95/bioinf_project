//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_IVITERBI_H
#define BIOINF_IVITERBI_H

typedef unsigned char byte;

#include <string>
#include "Sequence.h"

class IViterbi {

protected:
    const float (&emission_probabilities)[5][5];
    const float (&trans_prob)[3][3];

    std::map<char, int> _lookupTable;

public:
    // pure virtual function providing interface framework.
    virtual void alignSequences(Sequence *first, Sequence *second) = 0;

    IViterbi(const float (&emission_probabilities)[5][5], const float (&trans_prob)[3][3]) :
    emission_probabilities(emission_probabilities),
    trans_prob(trans_prob) {

        this->_lookupTable.insert(std::pair<char, int>('A', 0));
        this->_lookupTable.insert(std::pair<char, int>('C', 1));
        this->_lookupTable.insert(std::pair<char, int>('G', 2));
        this->_lookupTable.insert(std::pair<char, int>('T', 3));
        this->_lookupTable.insert(std::pair<char, int>('-', 4));


        this->_lookupTable.insert(std::pair<char, int>('d', 0));
        this->_lookupTable.insert(std::pair<char, int>('t', 1));
        this->_lookupTable.insert(std::pair<char, int>('e', 2));
    }

protected:
    float max(float v1, float v2, byte first, byte second, byte* result);
    float max(float m, float x, float y, byte* result);

    byte M = 1;
    byte X = 2;
    byte Y = 3;
};


#endif //BIOINF_IVITERBI_H
