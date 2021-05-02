import itertools
import copy
from warnings import warn


X, Y, Z, I = "X", "Y", "Z", "."


def dimensions(stab):

    width = len(stab[0][1])
    height = len(stab)

    return width, height


def display(stabilizer, shownumbers=False, logicals=None):

    width, height = dimensions(stabilizer)

    if shownumbers:
        print("Ph", end="\t")
        for k in range(width):
            print(k + 1, end="\t")

        print()

    print("=" * ((width + 1) * 4))

    for k, generator in enumerate(stabilizer):

        if logicals is not None:
            if k == logicals:
                print("- - " * (width + 1))

        if generator[0] == 0:
            print("+", end="\t")
        elif generator[0] == 1:
            print("+i", end="\t")
        elif generator[0] == 2:
            print("-", end="\t")
        elif generator[0] == 3:
            print("-i", end="\t")

        for w in generator[1]:
            print(w, end="\t")
        print()


def times(opx, opy):

    # if len(opx) > 2 or len(opy) > 2:

    #     if len(opx) == len(opy):

    #         length = len(opx)

    #         newphase = 0
    #         newout = [I] * length

    #         for k in range(length):

    #             kphase, kout = times((0, opx[k]), (0, opy[k]))
    #             newphase += kphase
    #             newout[k] = kout

    #         return tuple(newout.insert(0, newphase))

    #     else:
    #         raise ValueError("Can't times incompatible lengths", opx, opy)

    phx, x = opx
    phy, y = opy

    if x == I:
        return ((0 + phx + phy) % 4, y)
    if y == I:
        return ((0 + phx + phy) % 4, x)

    if (x, y) == (X, X) or (x, y) == (Y, Y) or (x, y) == (Z, Z):
        return ((0 + phx + phy) % 4, I)
    elif (x, y) == (X, Y):
        return ((1 + phx + phy) % 4, Z)
    elif (x, y) == (Y, X):
        return ((3 + phx + phy) % 4, Z)
    elif (x, y) == (Y, Z):
        return ((1 + phx + phy) % 4, X)
    elif (x, y) == (Z, Y):
        return ((3 + phx + phy) % 4, X)
    elif (x, y) == (Z, X):
        return ((1 + phx + phy) % 4, Y)
    elif (x, y) == (X, Z):
        return ((3 + phx + phy) % 4, Y)

    raise ValueError("Can't times values", opx, opy)


def genprod(gen1, gen2):

    phase = gen1[0] + gen2[0]
    paulis = []

    for w1, w2 in zip(gen1[1], gen2[1]):

        newphase, newpauli = times((0, w1), (0, w2))
        phase += newphase
        paulis.append(newpauli)

    return (phase % 4, paulis)


def commutes(gen1, gen2):

    sign = 0

    for w1, w2 in zip(gen1[1], gen2[1]):
        if w1 != I and w2 != I:
            if w1 != w2:
                sign += 1

    return sign % 2 == 0


def check(stab, logicals=0, verbose=False):

    assert logicals % 2 == 0

    for k1, logical in enumerate(stab[:logicals]):

        if k1 % 2 == 0:

            nextlogical = stab[k1 + 1]
            if commutes(logical, nextlogical):

                if verbose:
                    print("Logical pair", k1 + 1, "and", k1 + 2,
                          "commute\n")

                return False

        for k2, gen in enumerate(stab[logicals:]):
            if not commutes(logical, gen):

                if verbose:
                    print("Logical/generators number", k1 + 1, "and", k2 + 1,
                          "don't commute\n", gen1, "\n", gen2)

                return False

    for (k1, gen1), (k2, gen2) in itertools.combinations(enumerate(stab[logicals:]), 2):

        if not commutes(gen1, gen2):

            if verbose:
                print("Generators number", k1 + 1, "and", k2 + 1,
                      "don't commute\n", gen1, "\n", gen2)

            return False

    return True


def measure(stab, qubit, pauli, sign=+1, logicals=0, verbose=True):

    assert check(stab, verbose=True, logicals=logicals)
    assert sign == +1 or sign == -1

    if sign == +1:
        phase = 0
    elif sign == -1:
        phase = 2

    qid = qubit - 1
    newstab = copy.deepcopy(stab)
    width, height = dimensions(stab)

    measurer = (phase, [I] * width)
    measurer[1][qid] = pauli

    anticommerid = -1
    anticommer = None

    for k, gen in enumerate(newstab):

        if k > logicals - 1 and not commutes(gen, measurer):

            if anticommer is None:
                anticommer = gen
                anticommerid = k
                break

    if anticommer is not None:

        for k, gen in enumerate(newstab):

            if not commutes(gen, measurer):
                newstab[k] = genprod(gen, anticommer)

        newstab[anticommerid] = measurer

        for k, gen in enumerate(newstab):

            if gen[1][qid] == pauli and k != anticommerid:
                newstab[k] = genprod(gen, measurer)

    else:
        if verbose:
            print("Measurement commutes with stabilizer.")
        for k, gen in enumerate(newstab):

            if not commutes(gen, measurer):
                raise ValueError(
                    "Measurement commutes with stabilizer, but isn't in stabilizer")

            if gen[1][qid] == pauli:
                newstab[k] = genprod(gen, measurer)

    return newstab


def measurepattern(stab, pattern, verbose=False, logicals=0):

    newstab = copy.deepcopy(stab)

    if verbose:
        display(newstab, shownumbers=True, logicals=logicals)

    for meas in pattern:

        qubit, pauli, sign = meas

        if verbose:
            signstr = "+" if sign == +1 else "-"
            print(f"Measuring #{qubit} along {signstr}{pauli}")

        newstab = measure(newstab, qubit, pauli, sign, logicals)

        if verbose:
            display(newstab, logicals=logicals)

        assert check(newstab, logicals)

    return newstab


def apply_func(dict, stab, qubit1, qubit2):

    qid1, qid2 = qubit1 - 1, qubit2 - 1
    newstab = copy.deepcopy(stab)

    for k, gen in enumerate(newstab):

        phase, paulis = gen
        pauli1, pauli2 = paulis[qid1], paulis[qid2]

        xpart1 = dict[(X, I)] if pauli1 == X or pauli1 == Y else (0, (I, I))
        zpart1 = dict[(Z, I)] if pauli1 == Z or pauli1 == Y else (0, (I, I))
        xpart2 = dict[(I, X)] if pauli2 == X or pauli2 == Y else (0, (I, I))
        zpart2 = dict[(I, Z)] if pauli2 == Z or pauli2 == Y else (0, (I, I))

        newphase = 0

        if pauli1 == Y:
            newphase += 3
        if pauli2 == Y:
            newphase += 3

        addedgen = (0, (I, I))
        addedgen = genprod(addedgen, xpart1)
        addedgen = genprod(addedgen, zpart1)
        addedgen = genprod(addedgen, xpart2)
        addedgen = genprod(addedgen, zpart2)

        newphase += addedgen[0]
        newpauli1, newpauli2 = addedgen[1]

        newstab[k][0] += newphase
        newstab[k][1][qid1] = newpauli1
        newstab[k][1][qid2] = newpauli2

    return newstab


_CZ_dict = {
    (I, X): (0, (Z, X)),
    (I, Z): (0, (I, Z)),
    (X, I): (0, (X, Z)),
    (Z, I): (0, (Z, I))
}


def apply_cz(stab, qubit1, qubit2):
    return apply_func(_CZ_dict, stab, qubit1, qubit2)


def diag_stab(ops, padfront=0, padback=0):

    stab = []
    length = len(ops) + padfront + padback

    for k, op in enumerate(ops):
        paulis = [I] * length
        paulis[k + padfront] = op[1]
        stab.append([op[0], paulis])

    return stab


def entangle(stab, edges):

    newstab = copy.deepcopy(stab)

    for qubit1, qubit2 in edges:
        newstab = apply_cz(newstab, qubit1, qubit2)

    return newstab


def stab_with_inputs(inputs, qubits):

    msmts = qubits - inputs

    stab = diag_stab([(0, X)] * msmts, padfront=inputs)
    logicals = [[0, [I for _ in range(qubits)]] for _ in range(2 * inputs)]

    for k in range(inputs):

        logicals[2 * k][1][k] = X
        logicals[2 * k + 1][1][k] = Z

    stab = logicals + stab

    return stab


def entangled_square_with_inputs(inputs):

    return entangled_rectangle_with_inputs(inputs, 2 * inputs - 1)


def entangled_rectangle_with_inputs(inputs, width):

    height = 2 * inputs - 1
    qubits = 2 * inputs + height * width
    msmts = qubits - inputs

    stab = stab_with_inputs(inputs, qubits)

    topleft = inputs + 1
    bottomright = qubits - inputs

    # grid
    edges = [(i, i + 1) for firstone in range(topleft, bottomright, height)
             for i in range(firstone, firstone + height - 1)]
    edges += [(i, i + height)
              for i in range(topleft, bottomright - height + 1)]

    # entangle inputs and outputs
    edges += [(i, 2 * i + inputs - 1) for i in range(1, inputs + 1)]
    edges += [(bottomright - 2 * i, qubits - i) for i in range(inputs)]

    stab = entangle(stab, edges)

    return stab, msmts


# -----


def verify_cnot_rauss():

    INPUTS = 2
    QUBITS = 15

    stab = stab_with_inputs(INPUTS, QUBITS)

    edges = [(1, 3), (3, 5), (5, 7), (7, 10), (10, 12), (12, 14),
             (7, 8), (8, 9),
             (2, 4), (4, 6), (6, 9), (9, 11), (11, 13), (13, 15)]

    display(stab, logicals=4, shownumbers=True)

    print("Entangling")
    stab = entangle(stab, edges)

    pattern = [(1, X, +1),
               (2, X, +1),
               (3, Y, +1),
               (4, X, +1),
               (5, Y, +1),
               (6, X, +1),
               (7, Y, +1),
               (8, Y, +1),
               (9, Y, +1),
               (10, Y, +1),
               (11, X, +1),
               (12, Y, +1),
               (13, X, +1)]

    measurepattern(stab, pattern, verbose=True, logicals=4)


def verify_cnot_alt():

    INPUTS = 2
    QUBITS = 8

    stab = stab_with_inputs(INPUTS, QUBITS)

    edges = [(1, 5), (5, 6), (6, 7),
             (1, 3), (3, 4),
             (2, 4), (4, 8)]

    display(stab, logicals=4, shownumbers=True)

    print("Entangling")
    stab = entangle(stab, edges)

    pattern = [(1, X, -1),
               (2, X, +1),
               (3, Y, +1),
               (4, Y, +1),
               (5, Y, +1),
               (6, Y, +1), ]

    measurepattern(stab, pattern, verbose=True, logicals=4)


def verify_threelong_identity():

    INPUTS = 1
    QUBITS = 4

    stab = stab_with_inputs(INPUTS, QUBITS)

    edges = [(1, 2), (2, 3), (3, 4)]

    display(stab, logicals=2, shownumbers=True)

    print("Entangling")
    stab = entangle(stab, edges)

    pattern = [(1, Y, +1),
               (2, Y, +1),
               (3, Y, +1)]

    measurepattern(stab, pattern, verbose=True, logicals=2)


def verify_nswap(n):

    stab, msmts = entangled_square_with_inputs(n)

    display(stab, shownumbers=True, logicals=2 * n)

    pattern = [(i, X, +1) for i in range(1, msmts + 1)]

    measurepattern(stab, pattern, verbose=True, logicals=2 * n)


if __name__ == "__main__":

    pass
