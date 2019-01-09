#include <cstring>
#include <stdio.h>
#include <iostream>
#include "Viterbi.h"
#include "FastaParser.h"

typedef unsigned long ulong;

void Viterbi::alignSequences(Sequence *first, Sequence *second) {
    std::vector<char> first_sequence = first->getSequence();
    std::vector<char> second_sequence = second->getSequence();

    ulong n = first_sequence.size();
    ulong m = second_sequence.size();

    float* viterbi_match = new float[m];
    float* viterbi_insert_x = new float[m];
    float* viterbi_insert_y = new float[m];

    memset(viterbi_match, 0, sizeof(float) * m);
    memset(viterbi_insert_x, 0, sizeof(float) * m);
    memset(viterbi_insert_y, 0, sizeof(float) * m);

    viterbi_match[0] = 1;

    // ovdje zapisujemo prijelaze za sva stanja koja cemo koristiti u backtracku
    byte* transitions_m = new byte[n * m];
    byte* transitions_x = new byte[n * m];
    byte* transitions_y = new byte[n * m];

    memset(transitions_m, 0, sizeof(char) * n * m);
    memset(transitions_x, 0, sizeof(char) * n * m);
    memset(transitions_y, 0, sizeof(char) * n * m);

    byte previous;

    // spremanje medurezultata prethodnog koraka
    for (int i = 0; i < n; i++) {
        float* tmp_m = new float[m];
        float* tmp_x = new float[m];
        float* tmp_y = new float[m];

        memset(tmp_m, 0, sizeof(float) * m);
        memset(tmp_x, 0, sizeof(float) * m);
        memset(tmp_y, 0, sizeof(float) * m);

        float v1, v2, v3, max;

        const float v = trans_prob[0][0];

        char &c2 = second_sequence.at(0);
        char &c1 = first_sequence.at(i);

        int pos1 = lookup.at(c1);
        int pos2 = lookup.at(c2);

        float prob = emission_probabilities[pos1][pos2];

        tmp_m[0] = prob * v * viterbi_match[0];  // prijelaz iz M u M

        *(transitions_m + i * n + 0) = M;

        v1 = trans_prob[0][1] * viterbi_match[0]; // iz M u X
        v2 = trans_prob[1][1] * viterbi_insert_x[0]; // iz X u X
        max = this->max(v1, v2, M, X, &previous);

        tmp_x[0] = emission_probabilities[lookup.at(first_sequence.at(i))][lookup.at('-')] * max;
        *(transitions_x + i * n + 0) = previous;

        tmp_y[0] = 0; // za j = 1 ovo je uvijek nula prema pretpostavci
        *(transitions_y + i * n + 0) = M;

        for (ulong j = 1; j < m; j++) {
            // moramo vidjet sta cemo dalje zapisati u tmp polja, i sta cemo koristiti za racunanje

            // ovo sve uzimamo iz prethodne iteracije zbog (i - 1) u formuli
            v1 = trans_prob[0][0] * viterbi_match[j - 1]; // iz M u M
            v2 = trans_prob[1][0] * viterbi_insert_x[j - 1]; // iz X u M
            v3 = trans_prob[2][0] * viterbi_insert_y[j - 1]; // iz Y u M

            max = this->max(v1, v2, v3, &previous); // pronademo max

            tmp_m[j] = emission_probabilities[lookup.at(first_sequence.at(i))][lookup.at(second_sequence.at(j))] * max;
            *(transitions_m + i * n + j) = previous;

            v1 = trans_prob[0][1] * viterbi_match[j]; // iz M u X
            v2 = trans_prob[1][1] * viterbi_insert_x[j]; // iz X u X

            max = this->max(v1, v2, M, X, &previous);
            tmp_x[j] = emission_probabilities[lookup.at(first_sequence.at(i))][lookup.at('-')] * max;
            *(transitions_x + i * n + j)  = previous;

            v1 = trans_prob[0][2] * tmp_m[j - 1]; // iz M u Y
            v2 = trans_prob[2][2] * tmp_y[j - 1]; // iz Y u Y
            max = this->max(v1, v2, M, Y, &previous);

            tmp_y[j] = emission_probabilities[lookup.at('-')][lookup.at(second_sequence.at(j))] * max;
            *(transitions_y + i * n + j)  = previous;
        }

        memcpy(viterbi_match, tmp_m, m);
        memcpy(viterbi_insert_x, tmp_x, m);
        memcpy(viterbi_insert_y, tmp_y, m);

        delete [] tmp_m;
        delete [] tmp_x;
        delete [] tmp_y;
    }

    // svako slovo mapira se na prvi element 2D polja
    std::map<byte, byte*> transitions_lookup = {
            {M, transitions_m},
            {X, transitions_x},
            {Y, transitions_y}
    };

    // ovdje cemo zapisivat poravnanja
    std::vector<std::tuple<char, char>> aligned;

    // kad izademo iz petlji, u poljima na zadnjem mjestu ce biti vjerojatnoti svakog stanja
    // za vrijednosti n i m te od tuda zapocinjemo backtrack

    float vm = viterbi_match[m - 1];
    float vx = viterbi_insert_x[m - 1];
    float vy = viterbi_insert_y[m - 1];

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

    // ove dolje korake moramo ponavljat sve dok ne dobijemo ukupnu sekvencu stanja
    while (i > 0 && j > 0) {
        // dohvati vrijednost polja
        char state = *(trans + i * n + j);

        // iz nekog razloga mi switch neda da koristim M,X,Y
        // ovisno o stanju radimo odma poravnanje
        switch (state) {
            case 1: // M
                aligned.push_back(std::make_tuple(first_sequence.at(i), second_sequence.at(j)));

                j = j - 1;
                i = i - 1;
                break;
            case 2: // X --> emitira stanje xi,-
                aligned.push_back(std::make_tuple(first_sequence.at(i), '-'));
                i = i - 1;
                break;
            case 3: // Y --> emitira stanje -,yj
                aligned.push_back(std::make_tuple('-', second_sequence.at(j)));
                j = j - 1;
                break;
        }

        trans = transitions_lookup[state]; // iz onoga sta je zapisano u polju odredujemo sljedecu tablicu
    }

    delete [] transitions_x;
    delete [] transitions_y;
    delete [] transitions_m;

    delete [] viterbi_match;
    delete [] viterbi_insert_x;
    delete [] viterbi_insert_y;
}

Viterbi::Viterbi(const float (&emission_probabilities)[5][5], const float (&trans_prob)[3][3]): IViterbi(emission_probabilities, trans_prob) {}
