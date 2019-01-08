#include <cstring>
#include "Viterbi.h"

typedef unsigned long ulong;

float Viterbi::max(float v1, float v2) {
    if (v1 > v2) return v1;
    return v2;
}

float Viterbi::max(float v1, float v2, float v3) {
    float max = v1;
    if (v2 > max) max = v2;
    if (v3 > max) max = v3;
    return max;
}

std::string Viterbi::alignSequences(Sequence &first, Sequence &second) {
    std::vector<char> states;

    const std::vector<char> &first_sequence = first.getSequence();
    const std::vector<char> &second_sequence = second.getSequence();

    ulong n = first_sequence.size();
    ulong m = second_sequence.size();

    // trenutna pomocna polja
    float* viterbi_match = new float[m];
    float* viterbi_insert_x = new float[m];
    float* viterbi_insert_y = new float[m];

    memset(viterbi_insert_x, 0, sizeof(float) * m);
    memset(viterbi_insert_y, 0, sizeof(float) * m);
    memset(viterbi_match, 0, sizeof(float) * m);

    viterbi_match[0] = 1;

    const double tau = transition_probabilities[_lookupTable.at('t')];
    const double delta = transition_probabilities[_lookupTable.at('d')];
    const double epsilon = transition_probabilities[_lookupTable.at('e')];

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

        v1 = delta * viterbi_match[0];
        v2 = epsilon * viterbi_insert_x[0];
        max = this->max(v1, v2);

        tmp_x[0] = emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at('-')] * max;

        tmp_y[0] = 0; // za j = 1 ovo je uvijek nula prema pretpostavci

        for (ulong j = 1; j < m; j++) {

            // moramo vidjet sta cemo dalje zapisati u tmp polja, i sta cemo koristiti za racunanje

            // ovo sve uzimamo iz prethodne iteracije zbog (i - 1) u formuli
            v1 = (1 - 2 * delta - tau) * viterbi_match[j - 1];
            v2 = (1 - epsilon - tau) * viterbi_insert_x[j - 1];
            v3 = (1 - epsilon - tau) * viterbi_insert_y[j - 1];

            max = this->max(v1, v2, v3); // pronademo max

            tmp_m[j] = emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at(second_sequence.at(j))] * max;

            v1 = delta * viterbi_match[j];
            v2 = epsilon * viterbi_insert_x[j];

            max = this->max(v1, v2);
            tmp_x[j] = emission_probabilities[_lookupTable.at(first_sequence.at(i))][_lookupTable.at('-')] * max;


            v1 = delta * tmp_m[j - 1];
            v2 = epsilon * tmp_y[j - 1];
            max = this->max(v1, v2);

            tmp_y[j] = emission_probabilities[_lookupTable.at('-')][_lookupTable.at(second_sequence.at(j))] * max;
        }

        // izbrisi memoriju koju su zauzela stara polja
        delete [] viterbi_match;
        delete [] viterbi_insert_x;
        delete [] viterbi_insert_y;

        // zamijeni polja
        viterbi_match = tmp_m;
        viterbi_insert_x = tmp_x;
        viterbi_insert_y = tmp_y;

        tmp_m = NULL;
        tmp_x = NULL;
        tmp_y = NULL;
    }


    // Stari kod koji nije bas najbrzi
    /*for (unsigned long i = 0; i < n; i++) {
        for (unsigned long j = 0; j < m; j++) {

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
    }*/


    //TODO
    return "max probable path";


}

Viterbi::Viterbi(const double *transition_probabilities, const double **emission_probabilities) : IViterbi(
        transition_probabilities,
        emission_probabilities) {}

