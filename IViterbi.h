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
    const float (&emission_probabilities)[5][5];
    const float (&transition_probabilities)[3];
    std::map<char, int> lookup;
    const char gap = '-';
    const float tau;
    const float delta;
    const float epsilon;


public:
    // pure virtual function providing interface framework.
    virtual void alignSequences(Sequence *first, Sequence *second) = 0;

    IViterbi(const float (&transition_probabilities)[3], const float (&emission_probabilities)[5][5],
             std::map<char, int> &lookup) :
            emission_probabilities(emission_probabilities), transition_probabilities(transition_probabilities),
            lookup(lookup),
            tau(transition_probabilities[0]), delta(transition_probabilities[1]), epsilon(transition_probabilities[2]) {
    }

protected:
    double max(double v1, double v2, byte first, byte second, byte *result) {
        if (v1 > v2) {
            *result = first;
            return v1;
        }

        *result = second;
        return v2;
    };

    double max(double m, double x, double y, byte *result) {
        double max = this->max(m, x, M, X, result);

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
