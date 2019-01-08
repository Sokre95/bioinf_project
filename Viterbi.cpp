#include <cstring>
#include "Viterbi.h"

typedef unsigned long ulong;

float IViterbi::max(float v1, float v2, byte first, byte second, byte* result) {
    if (v1 > v2) {
        *result = first;
        return v1;
    }

    *result = second;
    return v2;
}

float IViterbi::max(float m, float x, float y, byte* result) {
    float max = this->max(m, x, M, X, result);

    if (max < y) {
        *result = Y;
        return y;
    }

    return max;
}

void IViterbi::alignSequences(Sequence &first, Sequence &second) {
    const std::vector<char> &first_sequence = first.getSequence();
    const std::vector<char> &second_sequence = second.getSequence();

    ulong n = first_sequence.size();
    ulong m = second_sequence.size();

    float* viterbi_match = new float[m];
    float* viterbi_insert_x = new float[m];
    float* viterbi_insert_y = new float[m];

    memset(viterbi_insert_x, 0, sizeof(float) * m);
    memset(viterbi_insert_y, 0, sizeof(float) * m);
    memset(viterbi_match, 0, sizeof(float) * m);

    viterbi_match[0] = 1;

    const float tau = transition_probabilities[_lookupTable.at('t')];
    const float delta = transition_probabilities[_lookupTable.at('d')];
    const float epsilon = transition_probabilities[_lookupTable.at('e')];

    // ovdje zapisujemo prijelaze za sva stanja koja cemo koristiti u backtracku
    byte transitions_m[n][m];
    byte transitions_x[n][m];
    byte transitions_y[n][m];

    memset(transitions_m, 0, sizeof(char) * n * m);
    memset(transitions_x, 0, sizeof(char) * n * m);
    memset(transitions_y, 0, sizeof(char) * n * m);

    byte previous;

    // spremanje medurezultata prethodnog koraka
    for (ulong i = 0; i < n; i++) {
        float* tmp_m = new float[m];
        float* tmp_x = new float[m];
        float* tmp_y = new float[m];

        float v1, v2, v3, max;

        memset(viterbi_insert_x, 0, sizeof(float) * m);
        memset(viterbi_insert_y, 0, sizeof(float) * m);
        memset(viterbi_match, 0, sizeof(float) * m);

        tmp_m[0] = emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at(
                second_sequence.at(0))] * (1 - 2*delta - tau) * viterbi_match[0];

        transitions_m[i][0] = M;

        v1 = delta * viterbi_match[0];
        v2 = epsilon * viterbi_insert_x[0];
        max = this->max(v1, v2, M, X, &previous);

        tmp_x[0] = emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at('-')] * max;
        transitions_x[i][0] = previous;

        tmp_y[0] = 0; // za j = 1 ovo je uvijek nula prema pretpostavci
        transitions_y[i][0] = M;

        for (ulong j = 1; j < m; j++) {
            // moramo vidjet sta cemo dalje zapisati u tmp polja, i sta cemo koristiti za racunanje

            // ovo sve uzimamo iz prethodne iteracije zbog (i - 1) u formuli
            v1 = (1 - 2 * delta - tau) * viterbi_match[j - 1];
            v2 = (1 - epsilon - tau) * viterbi_insert_x[j - 1];
            v3 = (1 - epsilon - tau) * viterbi_insert_y[j - 1];

            max = this->max(v1, v2, v3, &previous); // pronademo max

            tmp_m[j] = emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at(second_sequence.at(j))] * max;
            transitions_m[i][j] = previous;

            v1 = delta * viterbi_match[j];
            v2 = epsilon * viterbi_insert_x[j];

            max = this->max(v1, v2, M, X, &previous);
            tmp_x[j] = emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at('-')] * max;
            transitions_x[i][j] = previous;

            v1 = delta * tmp_m[j - 1];
            v2 = epsilon * tmp_y[j - 1];
            max = this->max(v1, v2, M, Y, &previous);

            tmp_y[j] = emission_probabilities[_lookupTable.at('-')][_lookupTable.at(second_sequence.at(j))] * max;
            transitions_y[i][j] = previous;
        }

        // izbrisi memoriju koju su zauzela stara polja
        delete [] viterbi_match;
        delete [] viterbi_insert_x;
        delete [] viterbi_insert_y;

        // zamijeni polja
        viterbi_match = tmp_m;
        viterbi_insert_x = tmp_x;
        viterbi_insert_y = tmp_y;
    }

    // svako slovo mapira se na prvi element 2D polja
    std::map<byte, byte*> transitions_lookup = {
            {M, transitions_m[0]},
            {X, transitions_x[0]},
            {Y, transitions_y[0]}
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
        trans = transitions_m[0];
    } else if (vx > vm && vx > vy) {
        trans = transitions_x[0];
    } else {
        trans = transitions_y[0];
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
}

Viterbi::Viterbi(const float *transition_probabilities, const float **emission_probabilities) : IViterbi(
        transition_probabilities,
        emission_probabilities) {}

