#ifndef BIOINF_SEQUENCE_H
#define BIOINF_SEQUENCE_H

#include <string>
#include <vector>

class Sequence {

    std::string description;

    std::vector<char> sequence;

public:
    /*
     * Returns the description of an sequence.
     */
    const std::string &getDescription() const;

    /*
     * Returns the sequence contents.
     */
    const std::vector<char> &getSequence() const;

    /*
     *Sets the sequence contents.
     */
    void setSequence(const std::vector<char> &sequence);

    /*
     *Sets the sequence description.
     */
    void setDescription(const std::string& description);
};


#endif //BIOINF_SEQUENCE_H
