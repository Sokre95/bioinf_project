//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_VITERBILOGODDS_H
#define BIOINF_VITERBILOGODDS_H


#include <map>
#include "IViterbi.h"

class ViterbiLogOdds : public IViterbi {

private:
    const double termination_constant_c;
    const double eta;

    const float (&transition_probabilities)[3];

    double s(char xi, char yj);

public:
    void alignSequences(Sequence *first, Sequence *second) override;

    ViterbiLogOdds(const float (&transition_probabilities)[3], const float (&emission_probabilities)[5][5],
                   std::map<char, int> &lookup,
                   double eta);


};


#endif //BIOINF_VITERBILOGODDS_H
