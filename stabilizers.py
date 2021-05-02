import numpy as np
from collections import namedtuple
import itertools

Stabilizer = namedtuple("Stabilizer", ["phases", "paulis"])


BINARY = {
    "I": 0,
    "Z": 1,
    "X": 2,
    "Y": 3,
    "+": 0,
    "+i": 1,
    "-": 2,
    "-i": 3
}


pcode = [["+", "X", "Z", "I"],
         ["+", "Z", "X", "Z"],
         ["+", "I", "Z", "X"]]


def powerset(iterable):
    """
    Stolen from the itertools doc page
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s) + 1))


def badd(x, y):
    """Does addition in the binary representation"""

    return np.bitwise_xor(x, y)


def bsum(iterable):

    if len(iterable) == 0:
        return 0
    else:
        x = iterable[0]
        for y in iterable[1:]:
            x = badd(x, y)
        return x


def commutes(brow1, brow2):
    """NB: takes unphased rows"""

    parity = True

    for i, (el1, el2) in enumerate(zip(brow1, brow2)):

        if el1 == 0 or el2 == 0 or el1 == el2:
            pass
        else:
            parity = not parity

    return parity


def stabilizer(bcode, checkcommutes=True, checklinindep=True):

    width = len(bcode[0]) - 1
    height = len(bcode)

    if height > width:
        raise ValueError("Stabilizer bcode can't be taller than wide")

    phases = np.zeros((height,), dtype=int)
    paulis = np.zeros((height, width), dtype=int)  # paulis[row, col]

    for i, row in enumerate(bcode):

        phases[i] = row[0]
        paulis[i] = row[1:]

    if checkcommutes:
        for brow1, brow2 in itertools.combinations(paulis, 2):
            if not commutes(brow1, brow2):
                raise ValueError("Stabilizer must have commuting elements")

    if checklinindep:

        for subset in powerset(paulis):
            if len(subset) > 0:
                if not np.any(bsum(subset)):
                    raise ValueError(
                        "Stabilizer must have linearly independent elements"
                    )

    return Stabilizer(phases, paulis)


def dimensions(stabilizer):

    return stabilizer.paulis.size


def binary_code(pcode):

    bcode = []

    for row in pcode:

        brow = [BINARY[row[0]]]

        for i, p in enumerate(row):
            if i > 0:
                brow.append(BINARY[p])

        bcode.append(brow)

    return np.array(bcode)
