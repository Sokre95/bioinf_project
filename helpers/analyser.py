import sys

from Bio.SubsMat import MatrixInfo

if len(sys.argv) != 3:
    print("You must provide files with aligned sequences!")
    exit(1)

file_1 = sys.argv[1]
file_2 = sys.argv[2]


def get_alignment(file):
    seq_1 = ""
    seq_2 = ""
    with open(file) as input_file:
        input_file.readline()
        first_sequence = input_file.readline().strip()
        input_file.readline()
        second_sequence = input_file.readline().strip()

        for j in range(len(first_sequence)):
            seq_1 += first_sequence[j]
            seq_2 += second_sequence[j]
    return seq_1, seq_2


def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]


def score_pairwise(seq1, seq2, matrix, gap_s, gap_e):
    score = 0
    gap = False
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gap_s
            else:
                score += score_match(pair, matrix)
        else:
            if '-' not in pair:
                gap = False
                score += score_match(pair, matrix)
            else:
                score += gap_e
    return score


seq_1_our, seq_2__our = get_alignment(file_1)
mafft_seq_1, mafft_seq_2 = get_alignment(file_2)

blosum = MatrixInfo.blosum62

print(score_pairwise(seq_1_our, seq_2__our, blosum, -5, -1))
print(score_pairwise(mafft_seq_1.upper(), mafft_seq_2.upper(), blosum, -5, -1))
