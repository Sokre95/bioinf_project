#ifndef BIOINF_SEQUENCE_H
#define BIOINF_SEQUENCE_H

#include <string>
#include <vector>

class Sequence {

    std::string description;

    std::vector<char> sequence;

public:

    const std::string &getDescription() const;

    const std::vector<char> &getSequence() const;

    void setSequence(const std::vector<char> &sequence);

    void setDescription(const std::string& description);
};


#endif //BIOINF_SEQUENCE_H
