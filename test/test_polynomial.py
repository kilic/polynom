from polynom.polynomial import evaluate, Polynomial, lagrange_interpolation
from polynom.ecc import Scalar


def test_cmp():
    A = Polynomial.rand(4)
    B = Polynomial.rand(3)
    assert A != B
    assert B != A
    assert A == A
    assert A != -A
    assert not A.is_zero()
    assert not B.is_zero()
    assert Polynomial.zero().is_zero()


def test_add():
    for i in range(1, 25):
        for j in range(1, 25):
            A = Polynomial.rand(i)
            B = Polynomial.rand(j)
            n = max(A.n(), B.n())
            Z = Polynomial.from_ints([0] * max(i, j))
            assert (A - A).is_zero()
            assert (A + B - A - B).is_zero()
            assert (A + B).n() == n
            assert A + Z == A
            assert Z + A == A
            assert Z + A + B == B + A + Z
            assert Z + A + B - A == Z + B
            assert A + B == A + B.coeffs
            assert A - B == A - B.coeffs
            assert B - A == -(A - B)
            assert B - A == -A + B
            C = Polynomial.rand(i)
            assert (A + B) + C == (C + B) + A
            u = Scalar.rand()
            assert A + u == A + Polynomial([u])
            assert A + u == A + [u]
            assert A - u == A - Polynomial([u])
            assert A - u == A - [u]


def test_eval():
    e = Scalar.rand()
    v = [Polynomial.rand(4) for _ in range(10)]
    vi0 = evaluate(e, *v)
    vi1 = [u(e) for u in v]
    for v0, v1 in zip(vi0, vi1):
        assert v0 == v1


def test_lagrange_interpolation():

    n = 10
    points = []
    for i in range(n):
        points.append((Scalar.rand(), Scalar.rand()))

    L_x = lagrange_interpolation(points)

    for i in range(n):
        assert (L_x(points[i][0]) == points[i][1])
