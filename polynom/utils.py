from polynom.scalar import Scalar, zero


def log2(n: int) -> int:
    assert n > 0
    r = 0
    while n > (1 << r):
        r += 1
    return r


def trim_zeros(coeffs: list[Scalar]) -> list[Scalar]:
    z = len(coeffs)
    for c in reversed(coeffs):
        if c != 0:
            break
        z -= 1
    return coeffs[:z]


def pad(u: list[Scalar], n: int, el=zero) -> list[Scalar]:
    return u[:] + (n - len(u)) * [el]


def bit_reverse(A: list[Scalar], n: int):
    l = len(A)
    B = A[:]
    assert l > 0
    assert ((l & (l - 1)) == 0)
    assert l == 1 << n
    for i in range(l):
        k = i
        r = 0
        for _ in range(n):
            r = (r << 1) | (k & 1)
            k >>= 1
        B[i] = A[r]
    return B


def i_to_fr(*v: int) -> list[Scalar]:
    return [e if isinstance(e, Scalar) else Scalar(e) for e in v]
