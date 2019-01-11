//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_IVITERBI_H
#define BIOINF_IVITERBI_H

typedef unsigned char byte;

#include <string>
#include <map>
#include "Sequence.h"

class IViterbi {

protected:
    const char gap = '-';
    const float tau;
    const float delta;
    const float epsilon;

    const float *transition_probabilities;
    float **emission_probabilities;
    std::map<char, int> lookup;


public:
    // pure virtual function providing interface framework.
    virtual void alignSequences(Sequence *first, Sequence *second) = 0;

    IViterbi(const float *transition_probabilities, float **emission_probabilities,
             std::map<char, int> &lookup) :
            emission_probabilities(emission_probabilities), transition_probabilities(transition_probabilities),
            lookup(lookup),
            delta(transition_probabilities[lookup['d']]), tau(transition_probabilities[lookup['t']]), epsilon(transition_probabilities[lookup['e']]) {
    }

protected:
    float max(float v1, float v2, byte first, byte second, byte *result) {
        if (v1 > v2) {
            *result = first;
            return v1;
        }

        *result = second;
        return v2;
    };

    float max(float m, float x, float y, byte *result) {
        float max = this->max(m, x, M, X, result);

        if (max < y) {
            *result = Y;
            return y;
        }

        return max;
    };

    byte M = 1;
    byte X = 2;
    byte Y = 3;
};


#endif //BIOINF_IVITERBI_H
