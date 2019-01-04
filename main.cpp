#include <iostream>
#include <fstream>

#include "FastaParser.h"

int main() {
    FastaParser parser("../files/sample.txt");

    const std::vector<Sequence *> sequences = parser.parse();

    for (int i = 0; i < sequences.size(); ++i) {
        Sequence *s = sequences.at(i);

        std::cout << s->getDescription() << std::endl;

        for (int j = 0; j < s->getSequence().size(); ++j) {
            std::cout << s->getSequence().at(j) << std::endl;
        }
    }

    return 0;
}
