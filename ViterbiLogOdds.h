//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_VITERBILOGODDS_H
#define BIOINF_VITERBILOGODDS_H


#include <map>
#include "IViterbi.h"

class ViterbiLogOdds : public IViterbi {
    double eta;

public:

    void alignSequences(Sequence *first, Sequence *second) override;

    ViterbiLogOdds(const float (&emission_probabilities)[5][5], const float (&trans_prob)[3][3], double eta);
};



#endif //BIOINF_VITERBILOGODDS_H
