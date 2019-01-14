#include "../include/ViterbiLogOddsOptimal.h"
#include <cmath>
#include <limits>
#include <cstring>
#include <iostream>

typedef unsigned long ulong;


float ViterbiLogOddsOptimal::s(char a, char b) {
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

void ViterbiLogOddsOptimal::alignSequences(Sequence *first, Sequence *second, std::vector<char>* top, std::vector<char>* bottom) {
    const std::vector<char> &first_sequence = first->getSequence();
    const std::vector<char> &second_sequence = second->getSequence();

    const float d = -std::log((delta * (1 - epsilon - tau)) / ((1 - eta) * (1 - 2 * delta - tau)));
    const float e = -std::log(epsilon / (1 - eta));

    ulong n = first_sequence.size();
    ulong m = second_sequence.size();

    auto *viterbi_match = new float[m];
    auto *viterbi_insert_x = new float[m];
    auto *viterbi_insert_y = new float[m];

    const float minus_infinity = std::numeric_limits<float>::lowest();

    auto **transitions_m = create2dArray<byte>(n, m);
    auto **transitions_x = create2dArray<byte>(n, m);
    auto **transitions_y = create2dArray<byte>(n, m);

    for (long i = 0; i < m; i++) {
        viterbi_match[i] = minus_infinity;
        viterbi_insert_x[i] = minus_infinity;
        viterbi_insert_y[i] = minus_infinity;
    }

    viterbi_match[0] = minus_infinity;
    viterbi_insert_x[0] = minus_infinity;

    byte result;
    float maxVal;

    for (ulong i = 0; i < n; i++) {
        if (this->printProgress && i % 14 == 0) {
            std::cout << "\r" << i << "/" << n << " lines processed" << std::flush;
        }

        auto *tmp_m = new float[m];
        auto *tmp_x = new float[m];
        auto *tmp_y = new float[m];

        if (i == 0) { // ovo je za slucaj i = 0, j = 0
            tmp_m[0] = -2 * std::log(eta);
            tmp_x[0] = minus_infinity;
            tmp_y[0] = minus_infinity;
        } else {
            tmp_m[0] = minus_infinity;
            transitions_m[i][0] = M; // mislim da je svejedno kaj ce tu bit jer nikad nece bit odabran

            // za sve druge i > 0, j = 0 procedura je sljedeca...

            // u tmp_x moramo izracunat kaj ce ic
            float vm = viterbi_match[0] - d;
            float vx = viterbi_insert_x[0] - e;

            maxVal = this->max(vm, vx, M, X, &result);
            transitions_x[i][0] = result;

            tmp_x[0] = maxVal;

            tmp_y[0] = minus_infinity;
            transitions_y[i][0] = M; // isto kao i gore, nikad nece bit odabran
        }

        for (ulong j = 1; j < m; j++) {
            float vm = viterbi_match[j - 1];
            float vx = viterbi_insert_x[j - 1];
            float vy = viterbi_insert_y[j - 1];

            maxVal = this->max(vm, vx, vy, &result);

            transitions_m[i][j] = result;

            float sVal = this->s(first_sequence.at(i), second_sequence.at(j));
            tmp_m[j] = sVal + maxVal;

            vm = viterbi_match[j] - d;
            vx = viterbi_insert_x[j] - e;

            maxVal = this->max(vm, vx, M, X, &result);

            transitions_x[i][j] = result;
            tmp_x[j] = maxVal;

            vm = tmp_m[j - 1] - d;
            vy = tmp_y[j - 1] - e;

            maxVal = this->max(vm, vy, M, Y, &result);

            transitions_y[i][j] = result;
            tmp_y[j] = maxVal;
        }

        delete[] viterbi_match;
        delete[] viterbi_insert_x;
        delete[] viterbi_insert_y;

        viterbi_match = tmp_m;
        viterbi_insert_x = tmp_x;
        viterbi_insert_y = tmp_y;
    }

    std::map<byte, byte **> transitions_lookup = {
            {M, transitions_m},
            {X, transitions_x},
            {Y, transitions_y}
    };

    double vm = viterbi_match[m - 1];
    double vx = viterbi_insert_x[m - 1];
    double vy = viterbi_insert_y[m - 1];

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


ViterbiLogOddsOptimal::ViterbiLogOddsOptimal(
        const float *transition_probabilities,
        float **emission_probabilities,
        std::map<char, int> &lookup, const float eta, bool printProgress)
        : IViterbi(transition_probabilities, emission_probabilities, lookup),
          eta(eta),
          termination_constant_c(std::log(1 - 2 * transition_probabilities[0] - transition_probabilities[1])
                                 - std::log(1 - transition_probabilities[2] - transition_probabilities[1])), printProgress(printProgress) {
}
