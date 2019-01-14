#include <iostream>
#include <fstream>
#include <chrono>
#include "FastaParser.h"
#include "Viterbi.h"
#include "ViterbiLogOdds.h"
#include "MleEstimator.h"

void print(long duration) {
    std::cout << std::endl;
    std::cout << "-----Kraj izvrsavanja programa-----" << std::endl;
    std::cout << "---------Analiza izvodenja---------" << std::endl;
    std::cout << "Vrijeme izvrsavanja: " << duration << " s" << std::endl;
}

int main() {

    FastaParser parser("../database/pairs/HIV1_REF_2010/p1.fasta");

    const std::vector<Sequence *> sequences = parser.parse();

    for (int i = 0; i < sequences.size(); ++i) {
        Sequence *s = sequences.at(i);

        std::cout << s->getDescription() << std::endl;

        for (int j = 0; j < s->getSequence().size(); ++j) {
            std::cout << s->getSequence().at(j) << std::flush;
        }
        std::cout << std::endl << std::endl;
    }

    std::cout << "Estimate:" << std::endl;
    auto *mleEstimator = new MleEstimator("../database/outputs_mafft/upcase/");
    mleEstimator->estimate();

    std::cout << "Viterbi:" << std::endl;
    auto *logOdds = new ViterbiLogOdds(mleEstimator->getAveragedTransitionProbabilities(),
                                                        mleEstimator->getEmissionProbabilities(),
                                                        mleEstimator->getLookupTable(), 0.01);

    std::vector<char> top;
    std::vector<char> bottom;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    logOdds->alignSequences(sequences.at(0), sequences.at(1), &top, &bottom);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    long duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

    // ispis...

    // ispis...

    /* for (auto it = top.rbegin(); it != top.rend(); ++it) {
         std::cout << *it << std::flush;
     }

     std::cout << std::endl;

     for (auto it = bottom.rbegin(); it != bottom.rend(); ++it) {
         std::cout << *it << std::flush;
     } */

    const int line_width = 80;
    ulong curr = top.size() - 1;

    bool isEnd;

    while (curr > 0) {
        int i;
        for (i = 0; i < line_width; i++) {
            if (curr == 0) {
                isEnd = true;
                break;
            }
            curr--;
            std::cout << top.at(curr) << std::flush;
        }

        std::cout << std::endl;

        if (isEnd) {
            curr += i;
        } else {
            curr += line_width - 1;
        }

        for (i = 0; i < line_width; i++) {
            if (curr == 0) break;
            curr--;
            std::cout << bottom.at(curr) << std::flush;
        }

        std::cout << std::endl;
        std::cout << std::endl;
    }

    print(duration);

    return 0;
}
