import sys
import math
from fractions import Fraction
from tqdm import tqdm
from random import *
from itertools import *

MAXN_EXHAUSTIVE = int(sys.argv[1]) # n <= 12 exhaustive counting

n=16
k=1
MASK = 2**n-1

def set_nk(n, k=1):
    globals()["n"] = n
    globals()["k"] = k
    globals()["MASK"] = 2**n - 1

def rol(x, r=1):
    r %= n
    x &= MASK
    y = (x << r) | (x >> (n - r))
    return y & MASK

def ror(x, r=1):
    r %= n
    return rol(x, n-r)

def add(x, y):
    # IMPORTANT TO DROP THE OUT CARRY
    return (x + y) & MASK

def rxdp_AL16(dx, dy, dz):
    # note: we return solution count to avoid with rationals
    # (factor 2^2n)
    dx >>= 1
    dy >>= 1
    dz >>= 1
    xor = dx ^ dy ^ dz
    neq = (dx ^ dy) | (dx ^ dz)

    maskn1 = MASK >> 1
    ISHLxor = (xor ^ (xor << 1)) & maskn1
    ISHLxor1 = (xor ^ (xor << 1) ^ 1) & maskn1
    SHLneq = (neq << 1) & maskn1
    d = bin(SHLneq).count("1")

    if ISHLxor1 & SHLneq == ISHLxor1:
        return 2**(2*n-d-3)
    if ISHLxor & SHLneq == ISHLxor:
        return 2**(2*n-d-3) * 3  #   3*2^-3 = 2^-1.415...
    return 0

def Tm(dx, dy, dz, w, m):
    # note: we return solution count to avoid with rationals
    # (factor 2^2m)
    xor = dx ^ dy ^ dz
    neq = (dx ^ dy) | (dx ^ dz)

    maskm = 2**m - 1
    SHLneq = (neq << 1) & maskm

    d = bin(SHLneq).count("1")
    res = 2**(2*m - d - 1)
    if xor in (0, 2**m-1):
        res += (-1)**w * 2**(2*m - m - 1)
    return res

# 2^(-d) + 2(-m)
# -----------------  = 1 +- 2^(d-m)
#       2(-d)

def rxdp_our(dx, dy, dz):
    # note: we return solution count to avoid with rationals
    # (factor 2^2n)
    xor = dx ^ dy ^ dz
    neq = (dx ^ dy) | (dx ^ dz)

    ISHLxor = (xor ^ (xor << 1)) & MASK
    SHLneq = (neq << 1) & MASK

    # validity criteria
    if ISHLxor & (SHLneq | 2**k | 1) != ISHLxor:
        return 0

    mask_k = (1 << k) - 1
    Tnk = Tm(dx >> k, dy >> k, dz >> k, xor & 1, n-k)
    Tk = Tm(dx & mask_k, dy & mask_k, dz & mask_k, (xor >> k) & 1, k)
    assert Tnk and Tk, (Tnk, Tk)
    return Tnk * Tk


def rxdp_exper(dx, dy, dz, progress=False):
    count_exper = 0
    xrange = tqdm(range(2**n)) if progress else range(2**n)
    for x in xrange:
        for y in range(2**n):
            left = rol((x+y) & MASK)
            right = ((rol(x) ^ dx) + (rol(y) ^ dy)) & MASK
            dz_test = left ^ right
            count_exper += dz == dz_test
    return count_exper


def example_1():
    print("Verifying Example 1. (Incorrect [AL16] value, correction 1/2)")
    set_nk(n=16, k=1)

    dx = 0b1001000000100101
    dy = 0b0010001010010010
    dz = 0b0100110101001000
    chi = 0b1111111111111111
    nu = 0b1111111111111111

    assert chi == dx ^ dy ^ dz
    assert nu == (dx ^ dy) | (dx ^ dz)

    mask_nu = (MASK >> 1) ^ 1  # drop 1 MSB and 1 LSB
    d = bin( nu & mask_nu ).count("1")
    assert d == n-2 == 14

    count_our = rxdp_our(dx, dy, dz)
    count_AL16 = rxdp_AL16(dx, dy, dz)
    assert count_our == 2**(32-18)
    assert count_AL16 == 2**(32-17)
    print("count our", count_our)
    print("count [AL16]", count_AL16)

    print("Exhaustive computation... this may take time (10 sec on pypy3)")
    count_exper = rxdp_exper(dx, dy, dz, progress=True)
    print("count experimental", count_exper)
    assert count_exper == count_our
    print("Done!\n")


def theorem_5_correction_factor():
    print("Random 1000 samples for each 3 <= n <= 32 (exhaustive check for n <= 12)")
    print("+ 1000 samples with forced dz")
    for n in range(3, 33):
        set_nk(n=n, k=1)
        print("n =", n)
        hits_agree = 0
        hits_disagree = 0
        for itr in tqdm(range(2000)):
            dx = randrange(2**n)
            dy = randrange(2**n)
            dz = randrange(2**n)
            if itr >= 1000:
                dz = dx ^ dy ^ (randrange(2) * MASK)
                dz ^= randrange(4)
                # note: may be invalid transition (count = 0)
                # so "hit disagree" won't trigger" >1000 times

            chi = dx ^ dy ^ dz
            nu = (dx ^ dy) | (dx ^ dz)

            mask_nu = (MASK >> 1) ^ 1  # drop 1 MSB and 1 LSB
            d = bin( nu & mask_nu ).count("1")

            count_our = rxdp_our(dx, dy, dz)
            count_AL16 = rxdp_AL16(dx, dy, dz)
            if n <= MAXN_EXHAUSTIVE:
                count_exper = rxdp_exper(dx, dy, dz, progress=False)
                assert count_our == count_exper

            # validity criteria agrees
            assert bool(count_AL16) == bool(count_our)
            if not count_AL16:
                continue

            # normal case
            if (chi >> 1) not in (0, MASK >> 1):
                assert count_our == count_AL16
                hits_agree += 1
            # corrected case
            else:
                # rephrased to avoid rationals
                to_add = count_AL16 * (-1)**(chi & 1) // 2**(n-1-d)
                assert count_our == count_AL16 + to_add
                hits_disagree += 1

        print("Hit cases: agree", hits_agree, "disagree", hits_disagree)
    print("Done!\n")


def sample_wrong_transition(sign):
    # randomly choose 1...1 or 0...0 for chi_L'
    chi = randrange(2) * (2**n-1)
    chi = 1 * (2**n-1)
    case = (chi >> 1) & 1 # 1...1 or 0...0

    # set chi_0 according to sign
    if chi & 1 != sign:
        chi ^= 1

    dx = randrange(2**n)
    dy = randrange(2**n)
    dz = (dx ^ dy ^ chi) & 1
    for i in range(1, n):
        # for case = 0  (case 1 symmetric)
        # 00 -> random of the 3
        if (dx >> i) & 1 == (dy >> i) & 1 == case:
            which = randrange(3)

            dx ^= 1 << i
            dy ^= 1 << i
            if case == 0:
                dz ^= 1 << i

            if which == 0:
                dx ^= 1 << i
            elif which == 1:
                dy ^= 1 << i
            elif which == 2:
                dz ^= 1 << i

        # for case = 0
        # 01 -> 011
        # 10 -> 101
        # 11 -> 110
        else:
            a = (dx >> i) & 1
            b = (dy >> i) & 1
            dz |= (a ^ b ^ case) << i
    return dx, dy, dz


def proposition_6():
    print("Random transitions of the shape given by Proposition 6 (correction cases)")
    print("1000 samples for each 3 <= n <= 32 (exhaustive check for n <= 12)")

    for n in range(3, 33):
        set_nk(n=n, k=1)
        print("n =", n)

        for itr in range(1000):
            sign = randrange(2)

            dx, dy, dz = sample_wrong_transition(sign)
            count_our = rxdp_our(dx, dy, dz)
            count_AL16 = rxdp_AL16(dx, dy, dz)
            if sign == 0:
                assert count_our == count_AL16 * 3 // 2
            if sign == 1:
                assert count_our == count_AL16 // 2

            if n <= MAXN_EXHAUSTIVE:
                count_exper = rxdp_exper(dx, dy, dz, progress=False)
                assert count_our == count_exper
    print("Done!\n")


def proposition_7():
    print("Counting valid transitions")

    prev = None
    for n in range(2, 9):
        for k in range(1, n):
            set_nk(n, k)
            print("n =", n, "k =", k)

            cnt = [0, 0]
            for dx in range(2**n):
                for dy in range(2**n):
                    for dz in range(2**n):
                        if rxdp_our(dx, dy, dz):
                            nu = (dx ^ dy) | (dx ^ dz)
                            msb_neq = nu >> (n-1)
                            cnt[msb_neq] += 1

            print("c0 c1", cnt, "sum", sum(cnt))
            assert sum(cnt) == 64 * 7**(n-2)
    print("Done!\n")


def corollary_2():
    print("(Corollary 2) Computing average probability (exhaustive)")

    # note: sampling dx, dy first uniformly, and then choosing dz
    # won't work, since different dx, dy have different numbers of valid trails
    # (even though they all have same probability contribution = 1)
    # and we compute an average probability over the set of all valid trails
    for n in range(2, 9):
        for k in range(1, n):
            set_nk(n, k)
            print("n =", n, "k =", k)
            total = 0
            n_valid = 0
            for dx in range(2**n):
                for dy in range(2**n):
                    cur = 0
                    for dz in range(2**n):
                        cnt = rxdp_our(dx, dy, dz)
                        if cnt:
                            total += cnt
                            n_valid += 1
                        cur += cnt
                    assert cur == 4**n

            prob = Fraction(total) / n_valid / 4**n  # count -> prob factor
            expec = Fraction(49, 64) * Fraction(4, 7)**n
            print("   ", "average", prob, "expected", expec)
            assert prob == expec
    print("Done!\n")

if __name__ == '__main__':
    set_nk(n=16, k=1)

    theorem_5_correction_factor()
    example_1()
    proposition_6()
    proposition_7()
    corollary_2()
