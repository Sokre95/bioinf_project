//
// Created by tomo on 1/8/19.
//

#include <cstring>
#include "Viterbi.h"

std::string Viterbi::alignSequences(Sequence &first, Sequence &second) {
    std::vector<char> first_sequence = first.getSequence();
    std::vector<char> second_sequence = second.getSequence();

    unsigned long n = first_sequence.size();
    unsigned long m = second_sequence.size();

    double viterbi_match[n][m];
    double viterbi_insert_x[n][m];
    double viterbi_insert_y[n][m];

    viterbi_match[0][0] = 1;
    memset(viterbi_insert_x, 0, sizeof(viterbi_insert_x[0][0]) * n * m);
    memset(viterbi_insert_y, 0, sizeof(viterbi_insert_y[0][0]) * n * m);

    std::vector<char> states;

    for (unsigned long i = 0; i < n; i++) {
        for (unsigned long j = 0; j < m; j++) {

            double tau = transition_probabilities[_lookupTable.at('t')];
            double delta = transition_probabilities[_lookupTable.at('d')];
            double epsilon = transition_probabilities[_lookupTable.at('e')];


            double temp = (1 - epsilon - tau);


            double max = std::max((1 - 2 * delta - tau) * viterbi_match[i - 1][j - 1],
                                  temp * viterbi_insert_x[i - 1][j - 1]);

            max = std::max(max, temp * viterbi_insert_y[i - 1][j - 1]);

            viterbi_match[i][j] = emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at(
                    second_sequence.at(j))] * max;


            viterbi_insert_x[i][j] =
                    emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at('-')] *
                    std::max(delta * viterbi_match[i - 1][j], epsilon * viterbi_insert_x[i - 1][j]);


            viterbi_insert_y[i][j] =
                    emission_probabilities[_lookupTable.at('-')][_lookupTable.at(second_sequence.at(j))] *
                    std::max(delta * viterbi_match[i][j - 1], epsilon * viterbi_insert_y[i][j - 1]);


        }
    }
    //TODO
    return "max probable path";


}

Viterbi::Viterbi(double *transition_probabilities, double **emission_probabilities) {
    this->emission_probabilities = emission_probabilities;
    this->transition_probabilities = transition_probabilities;

    this->_lookupTable.insert(std::pair<char, int>('A', 0));
    this->_lookupTable.insert(std::pair<char, int>('C', 1));
    this->_lookupTable.insert(std::pair<char, int>('G', 2));
    this->_lookupTable.insert(std::pair<char, int>('T', 3));
    this->_lookupTable.insert(std::pair<char, int>('-', 4));


    this->_lookupTable.insert(std::pair<char, int>('d', 0));
    this->_lookupTable.insert(std::pair<char, int>('t', 1));
    this->_lookupTable.insert(std::pair<char, int>('e', 2));

}

