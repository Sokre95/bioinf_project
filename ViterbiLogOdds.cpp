//
// Created by tomo on 1/8/19.
//

#include "ViterbiLogOdds.h"
#include <cmath>
#include <limits>
#include <cstring>
#include <iostream>

typedef unsigned long ulong;

double IViterbi::max(double v1, double v2, byte first, byte second, byte* result) {
    if (v1 > v2) {
        *result = first;
        return v1;
    }

    *result = second;
    return v2;
}

double IViterbi::max(double m, double x, double y, byte* result) {
    double max = this->max(m, x, M, X, result);

    if (max < y) {
        *result = Y;
        return y;
    }

    return max;
}

double ViterbiLogOdds::s(char a, char b) {
    double pab = emission_probabilities[lookup[a]][lookup[b]];
    double qa = emission_probabilities[lookup[a]][lookup[gap]];
    double qb = emission_probabilities[lookup[gap]][lookup[b]];

    double firstPartial = log(pab / (qa * qb));
    double secondPartial = log((1 - 2 * delta - tau) / ((1 - eta) * (1 - eta)));

    return firstPartial + secondPartial;
}

template< typename T> T arr_get(T* arr, ulong i, ulong j, ulong n) {
    return *(arr + i * n + j);
}

template<typename T> void arr_set(T* arr, ulong i, ulong j, ulong n, T value) {
    *(arr + i * n + j) = value;
}

void initialize(ulong n, ulong m, double* vm, double* vx, double* vy, double value) {
    for (ulong i = 1; i < n + 1; i++) {
        arr_set(vm, i, 0, m + 1, value);
        arr_set(vx, i, 0, m + 1, value);
        arr_set(vy, i, 0, m + 1, value);
    }

    for (ulong j = 1; j < m + 1; j++) {
        arr_set(vm, 0, j, m + 1, value);
        arr_set(vx, 0, j, m + 1, value);
        arr_set(vy, 0, j, m + 1, value);
    }
}

void ViterbiLogOdds::alignSequences(Sequence *first, Sequence *second) {
    const std::vector<char> &first_sequence = first->getSequence();
    const std::vector<char> &second_sequence = second->getSequence();

    double d = -log((delta * (1 - epsilon - tau)) / ((1 - eta) * (1 - 2 * delta - tau)));
    double e = -log(epsilon / (1 - eta));

    ulong n = first_sequence.size();
    ulong m = second_sequence.size();

    // pokazivaci na 2D polja
    auto* viterbi_match = new double[ (n + 1) * (m + 1) ]; //   [n + 1][m + 1];
    auto* viterbi_insert_x = new double[ (n + 1) * (m + 1) ]; //[n + 1][m + 1];
    auto* viterbi_insert_y = new double[ (n + 1) * (m + 1) ]; // [n + 1][m + 1];

    double minus_infinity = std::numeric_limits<double>::lowest();

    initialize(n, m, viterbi_match, viterbi_insert_x, viterbi_insert_y, minus_infinity);

    double val = -2 * log(eta);

    arr_set(viterbi_match, 1, 1, m + 1, val);

    double val1 = arr_get(viterbi_match, 1, 1, m+1);

    arr_set(viterbi_insert_x, 1, 1, m + 1, minus_infinity);
    arr_set(viterbi_insert_y, 1, 1, m + 1, minus_infinity);

    // inicijalizacija polja za pracenje stanja

    byte* transitions_m = new byte[n * m];
    byte* transitions_x = new byte[n * m];
    byte* transitions_y = new byte[n * m];

    // nisam siguran dal je ovo potrebno, al nek za sad ostane
    memset(transitions_m, 0, sizeof(char) * n * m);
    memset(transitions_x, 0, sizeof(char) * n * m);
    memset(transitions_y, 0, sizeof(char) * n * m);

    byte result;
    double maxVal;

    for (ulong i = 1; i < n + 1; i++) {
        for (ulong j = 1; j < m + 1; j++) {
            if (i == 1 && j == 1) continue;
            double vm = arr_get(viterbi_match, i - 1, j - 1, m+1);
            double vx = arr_get(viterbi_insert_x, i - 1, j - 1, m+1);
            double vy = arr_get(viterbi_insert_y, i - 1, j - 1, m+1);

            maxVal = this->max(vm, vx, vy, &result);

            arr_set(transitions_m, i - 1, j - 1, m, result);

            double sVal = this->s(first_sequence.at(i - 1), second_sequence.at(j - 1));
            arr_set(viterbi_match, i, j, m+1, (sVal + maxVal));

            vm = arr_get(viterbi_match, i - 1, j, m+1) - d;
            vx = arr_get(viterbi_insert_x, i - 1, j, m+1) - e;

            maxVal = this->max(vm, vx, M, X, &result);

            arr_set(transitions_x, i - 1, j - 1, m, result);
            arr_set(viterbi_insert_x, i, j, m+1, maxVal);


            vm = arr_get(viterbi_match, i, j - 1, m+1) - d;
            vy = arr_get(viterbi_insert_y, i, j - 1, m+1) - e;

            maxVal = this->max(vm, vy, M, Y, &result);

            arr_set(transitions_y, i - 1, j - 1, m, result);
            arr_set(viterbi_insert_y, i, j, m+1, maxVal);
        }
    }

    std::map<byte, byte*> transitions_lookup = {
            {M, transitions_m},
            {X, transitions_x},
            {Y, transitions_y}
    };

    double vm = arr_get(viterbi_match, n, m, m + 1);
    double vx = arr_get(viterbi_insert_x, n, m, m + 1);
    double vy = arr_get(viterbi_insert_y, n, m, m+ 1);

    // pokazivac na ispravno polje
    byte* trans;

    // indeksi po kojima idemo algoritmom backtrack
    ulong i = n - 1;
    ulong j = m - 1;

    if (vm > vx && vm > vy) {
        trans = transitions_m;
    } else if (vx > vm && vx > vy) {
        trans = transitions_x;
    } else {
        trans = transitions_y;
    }

    std::vector<std::tuple<char, char>> aligned;

    std::vector<char> top;
    std::vector<char> bottom;
    std::vector<char> states;

    // ove dolje korake moramo ponavljat sve dok ne dobijemo ukupnu sekvencu stanja
    while (i > 0 && j > 0) {
        // dohvati vrijednost polja
        char state = arr_get(trans, i, j, n);
        states.push_back(state);

        // iz nekog razloga mi switch neda da koristim M,X,Y
        // ovisno o stanju radimo odma poravnanje
        switch (state) {
            case 1: // M
                top.push_back(first_sequence.at(i));
                bottom.push_back(second_sequence.at(j));
                aligned.push_back(std::make_tuple(first_sequence.at(i), second_sequence.at(j)));

                j = j - 1;
                i = i - 1;
                break;
            case 2: // X --> emitira stanje xi,-
                top.push_back(first_sequence.at(i));
                bottom.push_back(gap);

                aligned.push_back(std::make_tuple(first_sequence.at(i), '-'));
                i = i - 1;
                break;
            case 3: // Y --> emitira stanje -,yj
                top.push_back(gap);
                bottom.push_back(second_sequence.at(j));

                aligned.push_back(std::make_tuple('-', second_sequence.at(j)));
                j = j - 1;
                break;
        }

        trans = transitions_lookup[state]; // iz onoga sta je zapisano u polju odredujemo sljedecu tablicu
    }

    ulong topL = top.size();
    ulong botL = bottom.size();

    for (ulong i = 0; i < top.size(); i++) {
        std::cout << top.at(i) << std::flush;
    }

    std::cout << std::endl;

    for (ulong i = 0; i < bottom.size(); i++) {
        std::cout << bottom.at(i) << std::flush;
    }


    delete [] viterbi_match;
    delete [] viterbi_insert_x;
    delete [] viterbi_insert_y;

    delete [] transitions_m;
    delete [] transitions_x;
    delete [] transitions_y;
}


ViterbiLogOdds::ViterbiLogOdds(
        const float (&transition_probabilities)[3],
        const float (&emission_probabilities)[5][5],
        const float (&trans_prob)[3][3],
        const double eta): IViterbi(emission_probabilities, trans_prob),
        transition_probabilities(transition_probabilities),
        eta(eta),
        delta(transition_probabilities[0]),
        tau(transition_probabilities[1]),
        epsilon(transition_probabilities[2]),
        termination_constant_c(log(1 - 2 * transition_probabilities[0] - transition_probabilities[1]) - log(1 - transition_probabilities[2] - transition_probabilities[1])) {
}