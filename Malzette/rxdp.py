from fractions import Fraction


def split(w, m):
    return w[:m], w[m:]


def compute_prob(k, dx, dy, dz):
    n = dx.n
    dxL, dxR = split(dx, n-k)
    dyL, dyR = split(dy, n-k)
    dzL, dzR = split(dz, n-k)
    num1 = compute(dxL, dyL, dzL, (dxR ^ dyR ^ dzR)[-1])
    num2 = compute(dxR, dyR, dzR, (dxL ^ dyL ^ dzL)[-1])
    return Fraction(num1 * num2, 4**n)


def compute(dx, dy, dz, lsb):
    xor = dx ^ dy ^ dz
    neq = (dx ^ dy) | (dx ^ dz)
    d = neq[1:].wt
    n = dx.n
    res = 2**(2*n-d-1)
    if 0 in (xor, ~xor):
        res += 2**(n-1) * (-1)**lsb
    return res
