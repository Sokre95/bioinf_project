import os
from collections import defaultdict

thymine = 0
guanine = 0
cytosine = 0
adenine = 0

number_of_transitions = 0
number_of_pairs = 0

pairs_frequency = defaultdict(int)
transition_frequency = defaultdict(int)
emission_frequency_x = defaultdict(int)
emission_frequency_y = defaultdict(int)
states = []

for filename in os.listdir("../outputs"):
    with open("../outputs/" + filename) as f:
        id_a = f.readline().strip('\n')
        sequence_a = f.readline().strip('\n')
        id_b = f.readline().strip('\n')
        sequence_b = f.readline().strip('\n')

        curr_states = ['B']
        for i in range(len(sequence_a)):
            first_seq = sequence_a[i]
            second_seq = sequence_b[i]

            number_of_pairs += 1
            if first_seq != '-' and second_seq != '-':
                curr_states.append('M')
                pairs_frequency[(first_seq, second_seq)] += 1
            if first_seq != '-' and second_seq == '-':
                curr_states.append('X')
                emission_frequency_x[(first_seq, second_seq)] += 1
            if second_seq != '-' and first_seq == '-':
                curr_states.append('Y')
                emission_frequency_y[(first_seq, second_seq)] += 1

        curr_states.append('E')
        states.append(curr_states)

for state in states:
    for i in range(len(state) - 1):
        number_of_transitions += 1
        transition_frequency[(state[i], state[i + 1])] += 1

print("Pairs probabilities: ")
for k, v in pairs_frequency.items():  # will become d.items() in py3k, 16.0
    print("p(%s) = %f\n" % (k, (v + 1.0) / (number_of_pairs + 2.0)))

print("Transitions probabilities: ")
for k, v in transition_frequency.items():  # will become d.items() in py3k, 14
    print("p(%s) = %f\n" % (k, (v + 1.0) / (number_of_transitions + 2.0)))

print("Emission probabilities x: ")
for k, v in emission_frequency_x.items():  # will become d.items() in py3k, 14
    print("p(%s) = %f\n" % (k, (v + 1.0) / (number_of_pairs + 2.0)))

print("Emission probabilities y: ")
for k, v in emission_frequency_y.items():  # will become d.items() in py3k, 14
    print("p(%s) = %f\n" % (k, (v + 1.0) / (number_of_pairs + 2.0)))
