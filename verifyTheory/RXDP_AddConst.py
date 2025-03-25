import time
from random import randint

#from tqdm import tqdm

# pip install binteger
from binteger import Bin  # Note: big endian


def LSB(x: Bin, i: int):
    """
    Returns the i-th LSB (from 0) of x.

    Helper to avoid the confusion of big-endian Bin.
    """
    return x[x.n-1-i]

# hardcoded transition matrices
Mrho = {
    (0,0,0): [2, 1, 0, 1],
    (0,0,1): [0, 0, 0, 0],
    (0,1,0): [1, 1, 0, 0],
    (0,1,1): [1, 0, 0, 1],
    (1,0,0): [2, 0, 0, 0],
    (1,0,1): [0, 1, 0, 1],
    (1,1,0): [1, 1, 0, 0],
    (1,1,1): [1, 0, 0, 1],
}

def carryProbAddConsts(delta, Delta, a, b):
    """
    (XDP over ModAdd with 2 constants)

    Counts C = #x : ((x ^ delta) + a)  ^  (x + b) = Delta,
    splitted into C0 + C1 where
    C1 counts only x's s.t. x + b >= 2^n

    Returns [C0, C1]
    """
    assert delta.n == Delta.n == a.n == b.n
    n = delta.n
    xor = delta ^ Delta ^ a ^ b
    T = 1, 0
    for i in range(1, n+1):
        va = LSB(a, i-1)
        vb = LSB(b, i-1)
        vd = LSB(delta, i-1)
        vD = LSB(Delta, i-1)
        if i < n:
            vxor_i = LSB(xor, i)
        else:
            vxor_i = 0
        rho = (
            vd ^ va ^ vb,
            vD ^ vd,
            vxor_i ^ va ^ vb,
        )
        if i < n:
            M = Mrho[rho]
        else:
            M = Mrho[0,0,0]  # [2 1 0 1]
        Ma,Mb,Mc,Md = M

        if vb:
            T = T[1], T[0]
        T = (
            Ma*T[0] + Mb*T[1],
            Mc*T[0] + Md*T[1],
        )
        if vb:
            T = T[1], T[0]
    return T


def rxdiffAddConst(k, delta, Delta, a):
    """
    RX-differential probability of addition with constant

    Counts C = #x : ((x0.rol(k) ^ delta) + a)  ^  (x0 + a).rol(k)  = Delta

    Returns C
    """
    n = delta.n
    assert delta.n == Delta.n == a.n
    assert 1 <= k <= n-1

    if 2*k > n:
        delta = delta.ror(k)
        Delta = Delta.ror(k)
        k = n- k

    aL = a[:k]
    aM = a[k:-k]
    aR = a[-k:]
    deltaL = delta[:-k]
    deltaR = delta[-k:]
    DeltaL = Delta[:-k]
    DeltaR = Delta[-k:]

    cR = LSB(deltaL ^ Bin.concat(aL, aM) ^ Bin.concat(aM, aR) ^ DeltaL, 0)
    cL = LSB(deltaR ^ aR ^ aL ^ DeltaR, 0)

    TL = carryProbAddConsts(deltaL, DeltaL, Bin.concat(aL, aM) + cR, Bin.concat(aM, aR))
    TR = carryProbAddConsts(deltaR, DeltaR, aL + cL, aR)
    return TL[cL] * TR[cR]


while True:

    for n in range(2, 17):
        print("n =", n)
        print("=======")

        print("test carryProbAddConsts")
        t0 = time.time()
        num1 = 0
        while time.time() < t0 + 1.0:
            x0 = Bin.random(n)
            delta = Bin.random(n)
            a = Bin.random(n)
            b = Bin.random(n)
            Delta = ((x0 ^ delta) + a) ^ (x0 + b)

            cnts = [0, 0]
            for x in Bin.iter(n):
                if ((x ^ delta) + a) ^ (x + b) == Delta:
                    cnts[x.int + b.int >= 2**n] += 1

            cnts2 = carryProbAddConsts(delta, Delta, a, b)
            assert cnts == list(cnts2)
            num1 += 1

        print("test rxdiffAddConst")
        t0 = time.time()
        num2 = 0
        while time.time() < t0 + 1.0:
            k = randint(1, n-1)
            x0 = Bin.random(n)
            delta = Bin.random(n)
            a = Bin.random(n)
            Delta = ((x0.rol(k) ^ delta) + a) ^ (x0 + a).rol(k)

            cnt = 0
            for x in Bin.iter(n):
                if ((x.rol(k) ^ delta) + a) ^ (x + a).rol(k) == Delta:
                    cnt += 1

            cnt2 = rxdiffAddConst(k, delta, Delta, a)
            assert cnt == cnt2
            num2 += 1

        print("did", num1, "+", num2, "tests")
        print()
