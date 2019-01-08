//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_VITERBI_H
#define BIOINF_VITERBI_H


#include <string>
#include <map>
#include "Sequence.h"
#include "IViterbi.h"

class Viterbi : public IViterbi {

public:

    void alignSequences(Sequence &first, Sequence &second) override;

    Viterbi(const float *transition_probabilities, const float **emission_probabilities);

private:
    float max(float v1, float v2, char first, char second, char* result);
    float max(float m, float x, float y, char* result);

    const char M = 1;
    const char X = 2;
    const char Y = 3;
};


#endif //BIOINF_VITERBI_H
