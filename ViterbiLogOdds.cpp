//
// Created by tomo on 1/8/19.
//

#include "ViterbiLogOdds.h"
#include <cmath>
#include <limits>

void ViterbiLogOdds::alignSequences(Sequence *first, Sequence *second) {
    const std::vector<char> &first_sequence = first->getSequence();
    const std::vector<char> &second_sequence = second->getSequence();

    unsigned long n = first_sequence.size();
    unsigned long m = second_sequence.size();

    double viterbi_match[n + 1][m + 1];
    double viterbi_insert_x[n + 1][m + 1];
    double viterbi_insert_y[n + 1][m + 1];

    //indexi pocinju od -1 pa je sve shiftano u desno za 1
    viterbi_match[1][1] = -2 * log(eta);
    viterbi_insert_x[1][1] = std::numeric_limits<double>::lowest();
    viterbi_insert_y[1][1] = std::numeric_limits<double>::lowest();

    for (int i = 1; i < n + 1; i++) {
        viterbi_insert_x[i][0] = std::numeric_limits<double>::lowest();
        viterbi_insert_y[i][0] = std::numeric_limits<double>::lowest();
        viterbi_match[i][0] = std::numeric_limits<double>::lowest();
    }

    for (int j = 1; j < m + 1; j++) {
        viterbi_insert_x[0][j] = std::numeric_limits<double>::lowest();
        viterbi_insert_y[0][j] = std::numeric_limits<double>::lowest();
        viterbi_match[0][j] = std::numeric_limits<double>::lowest();
    }


    const double tau = 0;// transition_probabilities[_lookupTable.at('t')];
    const double delta = 0;// transition_probabilities[_lookupTable.at('d')];
    const double epsilon = 0;// transition_probabilities[_lookupTable.at('e')];

    const double termination_constant_c = log(1 - 2 * delta - tau) - log(1 - epsilon - tau);


    for (unsigned long i = 1; i < n + 1; i++) {
        for (unsigned long j = 1; j < m + 1; j++) {
            if (i != 1 && j != 1) {


                double max = std::max(viterbi_match[i - 1][j - 1], viterbi_insert_x[i - 1][j - 1]);
                max = std::max(max, viterbi_insert_y[i - 1][j - 1]);

                double s = log(emission_probabilities[_lookupTable.at(first_sequence.at(i - 1))][_lookupTable.at(
                        second_sequence.at(j - 1))] /
                               (emission_probabilities[_lookupTable.at(first_sequence.at(i - 1))][_lookupTable.at(
                                       '-')] *
                                       emission_probabilities[_lookupTable.at('-')][_lookupTable.at(
                                        second_sequence.at(j - 1))])) +
                           log((1 - 2 * delta - tau) / ((1 - eta) * (1 - eta)));

                double d = -log((delta * (1 - epsilon - tau)) / ((1 - eta) * (1 - 2 * delta - tau)));
                double e = -log(epsilon / (1 - eta));

                viterbi_match[i][j] = s + max;
                viterbi_insert_x[i][j] = std::max(viterbi_match[i - 1][j] - d, viterbi_insert_x[i - 1][j] - e);
                viterbi_insert_y[i][j] = std::max(viterbi_match[i][j - 1] - d, viterbi_insert_y[i][j - 1] - e);
            }
        }
    }

}


ViterbiLogOdds::ViterbiLogOdds(const float (&emission_probabilities)[5][5], const float (&trans_prob)[3][3],
                               const double eta)
        : IViterbi(emission_probabilities, trans_prob), eta(eta) {}