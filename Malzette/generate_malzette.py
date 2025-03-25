import sys
import math
from random import sample, randrange, shuffle, seed
from fractions import Fraction

from binteger import Bin

from rxdp import compute_prob
from malzette_data import rotations

# reproducible in PyPy3
SEED = 2  # int(sys.argv[1])


def unselfxor(target, k, fix=True):
    """Invert the map target = c ^ rol(c, k)."""
    n = target.n
    g = math.gcd(k, n)

    err = 0
    for chunk in target.split(n=g):
        err ^= chunk.int
    err = Bin(err, n)
    if fix:
        target ^= err
    else:
        raise ValueError("No solution")

    res = target.list
    for i in range(g):
        for j in range(n//g):
            ii = (i+k) % n
            res[ii] ^= res[i]
            i = ii
    res = Bin(res).ror(k)

    assert res ^ res.rol(k) == target
    if fix:
        return res, err
    return res


def make_diffs(n, d, lsb=None):
    """Generate high-probability differentials
    with given `d` and possibly fixed lsbxor=`lsb`"""
    alpha = [0] * n
    beta = [0] * n
    delta = [0] * n

    poses = sample(list(range(n-1)), d) + [n-1]
    if lsb is None:
        xor = randrange(2)
    else:
        xor = lsb
    for i in range(n):
        if i not in poses:
            if xor is not None:
                alpha[i] = beta[i] = delta[i] = xor
            else:
                alpha[i] = beta[i] = delta[i] = xor = randrange(2)
        if i in poses:
            alpha[i] = randrange(2)
            beta[i] = randrange(2)
            delta[i] = randrange(2)
            if alpha[i] == beta[i] == delta[i]:
                delta[i] ^= 1
            if xor is not None and alpha[i] ^ beta[i] ^ delta[i] != xor:
                alpha[i] ^= 1
            if alpha[i] == beta[i] == delta[i]:
                beta[i] ^= 1
                delta[i] ^= 1
            xor = None

    dtest = 0
    for i in range(n-1):
        if alpha[i] == beta[i] == delta[i]:
            assert alpha[i+1] ^ beta[i+1] ^ delta[i+1] == alpha[i], i
        else:
            dtest += 1
    assert dtest == d
    alpha = Bin(alpha[::-1])
    beta = Bin(beta[::-1])
    delta = Bin(delta[::-1])
    return alpha, beta, delta


def generate_transitions(n, k, d, balanced_percentile=1/3, nbest=2):
    best = set()
    for _ in range(5000):
        aR, bR, cR = make_diffs(k, 0)
        aL, bL, cL = make_diffs(n-k, d)
        a = Bin.concat(aL, aR)
        b = Bin.concat(bL, bR)
        c = Bin.concat(cL, cR)
        pr = compute_prob(k, a, b, c)
        bal = abs((a.str + b.str + c.str).count("1") - 3*n//2)
        best.add((
            pr, bal,
            (a, b, c)
        ))

    bals = sorted({bal for _, bal, _ in best})
    bal_threshold = bals[int(len(bals) * balanced_percentile)]
    best = [(u, v, w) for u, v, w in best if v <= bal_threshold]

    probs = sorted({prob for prob, _, _ in best})
    pr_threshold = probs[-nbest] if nbest <= len(probs) else probs[0]
    print("probs", *["%.3f" % math.log(pr, 2) for pr in probs[-nbest:]])
    print("bals", sorted({bal for _, bal, _ in best}))

    best = [(u, v, w) for u, v, w in best if u >= pr_threshold - 0.0001]
    shuffle(best)
    print("samples", len(best), "unique")
    print(*best[0][2])
    print(*best[1][2])
    print(*best[2][2])
    print()
    return best


def inner_to_io(rno, inner):
    """Convert RX-differential over addition
    to Alzette's input/output differences"""
    a, b, c = inner
    xl = a
    xr = b.rol(rotations[rno][0])
    yl = c
    yr = xr ^ yl.ror(rotations[rno][1])
    return (xl, xr, yl, yr)


def maxrun(s):
    """Compute maximum run of a character in a string"""
    prev = None
    cur = 0
    mx = 1
    for c in s:
        if c == prev:
            cur += 1
            mx = max(mx, cur)
        else:
            prev = c
            cur = 1
    return mx


def main():
    n = 32
    k = 3

    rounds = 12

    from random import shuffle, seed, randrange

    for MALZETTE_VER in (1, 2):
        seed(0)
        diffs = [
            generate_transitions(n, k, d=0, balanced_percentile=1/3, nbest=3),
            generate_transitions(n, k, d=1, balanced_percentile=1/3, nbest=3),
            generate_transitions(n, k, d=2, balanced_percentile=1/3, nbest=3),
        ]

        if MALZETTE_VER == 1:
            # seed
            sd = SEED
            # d=0 to get best probability
            diffses = [diffs[0] for i in range(rounds)]
            quality = 0
        elif MALZETTE_VER == 2:
            # seed
            sd = SEED
            # alternate d=1 and d=0
            diffses = [diffs[(i+1)%2] for i in range(rounds)]
            quality = 1

        seed(sd)
        print("Malzette%d" % MALZETTE_VER)
        print("=========")
        print("seed", sd)

        rounds = 12

        rno = 0
        for dat in diffses[0]:
            _, _, diff = dat
            io = inner_to_io(rno=rno, inner=diff)
            seq = recurse(diffses, rounds, k, rno=1, dx=io[2], dy=io[3], seen=set(), quality=quality)
            if seq is not None:
                seq = [dat] + seq
                break

        consts = []
        rno = 0
        prob, _, diff = seq[0]
        print("Round", rno, "prob", math.log(prob, 2))
        print("diffs", *seq[0][2])
        total_prob = prob
        for rno in range(1, rounds):
            io_prev = inner_to_io(rno=rno-1, inner=seq[rno-1][2])
            io = inner_to_io(rno=rno, inner=seq[rno][2])
            dl = io_prev[2] ^ io[0]
            dr = io_prev[3] ^ io[1]
            cl, el = unselfxor(dl, k)
            cr, er = unselfxor(dr, k)
            assert el == er == 0
            # assert maxrun(cl.hex) <= 2 and maxrun(cr.hex) <= 2
            print("const", cl.hex, cr.hex, cl, cr)
            #print("     ", (cl^cl.rol(k)).hex, (cr^cr.rol(k)).hex)
            print()
            prob, _, diff = seq[rno]
            print("Round", rno, "prob", math.log(prob, 2))
            print("diffs", *diff)
            total_prob *= prob
            consts.append((cl, cr))

        # add random constant in the last round to keep the structure
        consts.append((Bin.random(32), Bin.random(32)))

        # compute the full differentials
        diff_in = inner_to_io(rno=0, inner=seq[0][2])[:2]
        diff_out = inner_to_io(rno=rounds-1, inner=seq[rounds-1][2])[2:]
        diff_out = (
            diff_out[0] ^ consts[-1][0] ^ consts[-1][0].rol(k),
            diff_out[1] ^ consts[-1][1] ^ consts[-1][1].rol(k),
        )
        print("const", consts[-1][0].hex, consts[-1][1].hex)
        print()

        print("RX-differential: k =", k)
        print(diff_in[0].hex, diff_in[1].hex)
        print(" --%dR-->" % rounds)
        print(diff_out[0].hex, diff_out[1].hex)
        print("total prob", math.log(total_prob, 2))
        print()
        print()

        if MALZETTE_VER == 1:
            consts0 = consts[:]
            diffs0 = seq[:]
            diff_in0 = diff_in[:]
            diff_out0 = diff_out[:]
        else:
            consts1 = consts[:]
            diffs1 = seq[:]
            diff_in1 = diff_in[:]
            diff_out1 = diff_out[:]

    with open("data.tex", "wt") as f:
        pr = lambda *args: print(*args, file=f)

        pr("Table constants\n\n")
        probs = [Fraction(1), Fraction(1)]
        for rno in range(rounds):
            pr(
                rno+1,
                "&",

                r"\texttt{%s}:\texttt{%s}" % (consts0[rno][0].hex, consts0[rno][1].hex),
                "&",

                "%.2f" % math.log(diffs0[rno][0], 2),
                "&",

                r"\texttt{%s}:\texttt{%s}" % (consts1[rno][0].hex, consts1[rno][1].hex),
                "&",

                "%.2f" % math.log(diffs1[rno][0], 2),
                #"&",

                r"\\"
            )
            probs[0] *= diffs0[rno][0]
            probs[1] *= diffs1[rno][0]

        pr(r"\bottomrule")
        pr(r"Total", "&", "&", "%.2f" % math.log(probs[0], 2), "&", "&", "%.2f" % math.log(probs[1], 2), r"\\")
        pr("\n\n")

        pr("Table trails\n\n")
        probs = [Fraction(1), Fraction(1)]
        for rno in range(rounds):
            pr(
                rno+1,
                "&",

                r"\fontsize{8}{8}\texttt{%s}:\texttt{%s}:\texttt{%s}" % (diffs0[rno][2][0].hex, diffs0[rno][2][1].hex, diffs0[rno][2][2].hex),
                "&",

                "%.2f" % math.log(diffs0[rno][0], 2),
                "&",

                r"\fontsize{8}{8}\texttt{%s}:\texttt{%s}:\texttt{%s}" % (diffs1[rno][2][0].hex, diffs1[rno][2][1].hex, diffs1[rno][2][2].hex),
                "&",

                "%.2f" % math.log(diffs1[rno][0], 2),
                #"&",

                r"\\"
            )
            probs[0] *= diffs0[rno][0]
            probs[1] *= diffs1[rno][0]
        pr(r"\bottomrule")
        pr(r"Total", "&", "&", "%.2f" % math.log(probs[0], 2), "&", "&", "%.2f" % math.log(probs[1], 2), r"\\")
        pr("\n\n")

        pr("consts_Malzette1 = [")
        for cl, cr in consts0:
            pr("    (0x%08x, 0x%08x)," % (cl, cr))
        pr("]")
        pr()
        pr("consts_Malzette2 = [")
        for cl, cr in consts1:
            pr("    (0x%08x, 0x%08x)," % (cl, cr))
        pr("]")
        pr()

        pr("differential_Malzette1 = [")
        for cl, cr in (diff_in0, diff_out0):
            pr("    (0x%08x, 0x%08x)," % (cl, cr))
        pr("]")
        pr("differential_Malzette2 = [")
        for cl, cr in (diff_in1, diff_out1):
            pr("    (0x%08x, 0x%08x)," % (cl, cr))
        pr("]")
        pr()

        pr("const uint32_t differential_trail_Malzette1[] = {")
        for _, _, (a, b, c) in diffs0:
            pr("    0x%08x," % c)
        pr("};")
        pr("const uint32_t differential_trail_Malzette2[] = {")
        for _, _, (a, b, c) in diffs1:
            pr("    0x%08x," % c)
        pr("};")
        pr()

        for cl, cr in (diff_in0, diff_out0, diff_in1, diff_out1):
            pr("%08x:%08x" % (cl, cr))


def recurse(diffses, rounds, k, rno, dx, dy, seen, quality=0):
    """Recursively search for a high-prob. trail using
    given round transition sets and constraints on constants
    quality=0 - no constraints on constants (Malzette1)
    quality=1 - strict constraints on constants (Malzette2)
    """
    if rno == rounds:
        return []
    lst = diffses[rno][:]
    shuffle(lst)
    for dat in lst:
        (_, _, diff) = dat
        io = inner_to_io(rno=rno, inner=diff)
        dl = dx ^ io[0]
        dr = dy ^ io[1]
        cl, el = unselfxor(dl, k)
        cr, er = unselfxor(dr, k)
        if el or er:
            continue

        if quality:
            if maxrun(cl.hex) >= 3 or maxrun(cr.hex) >= 3:
                continue

            if "0" in cr.hex:
                continue
            if "0" in cl.hex:
                continue
            if cl in seen or cr in seen:
                continue

        seen2 = seen.copy() | {cl, cr}
        sub = recurse(diffses, rounds, k, rno+1, io[2], io[3], seen2, quality=quality)
        if sub is not None:
            return [dat] + sub




if __name__ == '__main__':
    main()
