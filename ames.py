"""Enumerate two-block samples of Ames, Iowa. Python standard library only.
Inputs: Horvitz & Thompson (1952), Table 2, columns 2 and 3.
The sequential PPS design below uses the raw size proxy on each draw.
"""
from itertools import combinations
from math import isclose, sqrt
from statistics import mean

households = [19, 9, 17, 14, 21, 22, 27, 35, 20, 15, 18, 37, 12, 47, 27, 25, 25, 13, 19, 12]
sizes = [18, 9, 14, 12, 24, 25, 23, 24, 17, 14, 18, 40, 12, 30, 27, 26, 21, 9, 19, 12]
pairs = list(combinations(range(len(sizes)), 2))
total_size = sum(sizes)

def pps_pair(i, j):
    first_i = sizes[i] / total_size
    first_j = sizes[j] / total_size
    return (first_i * sizes[j] / (total_size - sizes[i])
            + first_j * sizes[i] / (total_size - sizes[j]))

def evaluate(label, probabilities, inverse):
    inclusion = [sum(p for (i, j), p in zip(pairs, probabilities) if k in (i, j))
                 for k in range(len(sizes))]
    assert isclose(sum(probabilities), 1) and isclose(sum(inclusion), 2)
    estimates = [sum(households[k] / inclusion[k] for k in pair) if inverse
                 else len(sizes) * mean(households[k] for k in pair) for pair in pairs]
    expected = sum(p * est for p, est in zip(probabilities, estimates))
    variance = sum(p * (est - expected)**2 for p, est in zip(probabilities, estimates))
    if inverse:
        assert isclose(expected, sum(households))
    print(f'{label}: expectation={expected:.3f}, bias={expected-sum(households):.3f}, SD={sqrt(variance):.3f}')

if __name__ == '__main__':
    equal = [1 / len(pairs)] * len(pairs)
    pps = [pps_pair(i, j) for i, j in pairs]
    evaluate('SRS, expanded mean', equal, False)
    evaluate('PPS, expanded mean', pps, False)
    evaluate('PPS, HT', pps, True)
