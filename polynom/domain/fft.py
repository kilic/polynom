from polynom.scalar import Scalar
from polynom.utils import bit_reverse, log2


def perform_fft(A: list[Scalar], domain: list[Scalar]) -> list[Scalar]:
    n = len(A)
    exp = log2(n)
    assert n == len(domain)
    assert n == 1 << exp
    A = bit_reverse(A, exp)
    d = n >> 1
    for s in range(1, exp + 1):
        m = 1 << s
        k = 0
        while k < n:
            mm = m >> 1
            for j in range(mm):
                w = domain[j * d]
                t = w * A[k + j + mm]
                u = A[k + j]
                A[k + j] = u + t
                A[k + j + mm] = u - t
            k += m
        d >>= 1
    return A