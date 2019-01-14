//
// Created by tomo on 1/8/19.
//

#include "ViterbiLogOdds.h"
#include <cmath>
#include <limits>
#include <cstring>
#include <iostream>

typedef unsigned long ulong;


float ViterbiLogOdds::s(char a, char b) {
    float pab = emission_probabilities[lookup[a]][lookup[b]];
    float qa = emission_probabilities[lookup[a]][lookup[gap]];
    float qb = emission_probabilities[lookup[gap]][lookup[b]];

    float firstPartial = std::log(pab / (qa * qb));
    float secondPartial = std::log((1 - 2 * delta - tau) / ((1 - eta) * (1 - eta)));

    return firstPartial + secondPartial;
}

template<typename T>
T **create2dArray(ulong n, ulong m) {
    auto **array = new T *[n];
    for (int i = 0; i < n; ++i)
        array[i] = new T[m];

    return array;
}


void ViterbiLogOdds::alignSequences(Sequence *first, Sequence *second, std::vector<char>* top, std::vector<char>* bottom) {
    const std::vector<char> &first_sequence = first->getSequence();
    const std::vector<char> &second_sequence = second->getSequence();

    const float d = -std::log((delta * (1 - epsilon - tau)) / ((1 - eta) * (1 - 2 * delta - tau)));
    const float e = -std::log(epsilon / (1 - eta));

    ulong n = first_sequence.size();
    ulong m = second_sequence.size();

    auto **viterbi_match = create2dArray<float>(n + 1, m + 1);
    auto **viterbi_insert_x = create2dArray<float>(n + 1, m + 1);
    auto **viterbi_insert_y = create2dArray<float>(n + 1, m + 1);

    float minus_infinity = std::numeric_limits<float>::lowest();

    for (ulong i = 1; i < n + 1; i++) {
        viterbi_match[i][0] = minus_infinity;
        viterbi_insert_x[i][0] = minus_infinity;
        viterbi_insert_y[i][0] = minus_infinity;
    }

    for (ulong j = 1; j < m + 1; j++) {
        viterbi_match[0][j] = minus_infinity;
        viterbi_insert_x[0][j] = minus_infinity;
        viterbi_insert_y[0][j] = minus_infinity;
    }

    float val = -2 * std::log(eta);

    viterbi_match[1][1] = val;

    viterbi_insert_x[1][1] = minus_infinity;
    viterbi_insert_y[1][1] = minus_infinity;

    auto **transitions_m = create2dArray<byte>(n, m);
    auto **transitions_x = create2dArray<byte>(n, m);
    auto **transitions_y = create2dArray<byte>(n, m);

    byte result;
    float maxVal;

    for (ulong i = 1; i < n + 1; i++) {
        for (ulong j = 1; j < m + 1; j++) {
            if (i == 1 && j == 1) continue;
            float vm = viterbi_match[i - 1][j - 1];
            float vx = viterbi_insert_x[i - 1][j - 1];
            float vy = viterbi_insert_y[i - 1][j - 1];

            maxVal = this->max(vm, vx, vy, &result);

            transitions_m[i - 1][j - 1] = result; //arr_set(transitions_m, i - 1, j - 1, m, result);

            float sVal = this->s(first_sequence.at(i - 1), second_sequence.at(j - 1));
            viterbi_match[i][j] = sVal + maxVal;

            vm = viterbi_match[i - 1][j] - d;
            vx = viterbi_insert_x[i - 1][j] - e;

            maxVal = this->max(vm, vx, M, X, &result);

            transitions_x[i - 1][j - 1] = result;
            viterbi_insert_x[i][j] = maxVal;


            vm = viterbi_match[i][j - 1] - d;
            vy = viterbi_insert_y[i][j - 1] - e;

            maxVal = this->max(vm, vy, M, Y, &result);

            transitions_y[i - 1][j - 1] = result;
            viterbi_insert_y[i][j] = maxVal;
        }
    }

    std::map<byte, byte **> transitions_lookup = {
            {M, transitions_m},
            {X, transitions_x},
            {Y, transitions_y}
    };

    double vm = viterbi_match[n][m];
    double vx = viterbi_insert_x[n][m];
    double vy = viterbi_insert_y[n][m];

    // pokazivac na ispravno polje
    byte **trans;

    // indeksi po kojima idemo algoritmom backtrack
    ulong i = n - 1;
    ulong j = m - 1;

    vx += termination_constant_c;
    vy += termination_constant_c;

    if (vm > vx && vm > vy) {
        trans = transitions_m;
    } else if (vx > vm && vx > vy) {
        trans = transitions_x;
    } else {
        trans = transitions_y;
    }

    // ove dolje korake moramo ponavljat sve dok ne dobijemo ukupnu sekvencu stanja
    while (i > 0 && j > 0) {
        // dohvati vrijednost polja
        char state = trans[i][j]; //arr_get(trans, i, j, n);

        // iz nekog razloga mi switch neda da koristim M,X,Y
        // ovisno o stanju radimo odma poravnanje
        switch (state) {
            case 1: // M
                top->push_back(first_sequence.at(i));
                bottom->push_back(second_sequence.at(j));
                j = j - 1;
                i = i - 1;
                break;
            case 2: // X --> emitira stanje xi,-
                top->push_back(first_sequence.at(i));
                bottom->push_back(gap);
                i = i - 1;
                break;
            case 3: // Y --> emitira stanje -,yj
                top->push_back(gap);
                bottom->push_back(second_sequence.at(j));
                j = j - 1;
                break;
        }

        trans = transitions_lookup[state]; // iz onoga sta je zapisano u polju odredujemo sljedecu tablicu
    }

    top->push_back(first_sequence.at(0));
    bottom->push_back(second_sequence.at(0));

    delete[] viterbi_match;
    delete[] viterbi_insert_x;
    delete[] viterbi_insert_y;

    delete[] transitions_m;
    delete[] transitions_x;
    delete[] transitions_y;
}


ViterbiLogOdds::ViterbiLogOdds(
        const float *transition_probabilities,
        float **emission_probabilities,
        std::map<char, int> &lookup, const float eta)
        : IViterbi(transition_probabilities, emission_probabilities, lookup),
          eta(eta),
          termination_constant_c(std::log(1 - 2 * transition_probabilities[0] - transition_probabilities[1])
                                 - std::log(1 - transition_probabilities[2] - transition_probabilities[1])) {
}