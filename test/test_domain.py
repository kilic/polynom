from polynom.scalar import Scalar, one
from polynom.polynomial import Polynomial
from polynom.domain.domain import Domain
from polynom.ecc.bn254.domain import domain_config, kx


def test_interpolation():
    for i in range(1, 5):
        domain = Domain(domain_config(i))

        A0 = Polynomial.rand(domain.n)
        A_x = domain.interpolate(A0)
        A1 = domain.evaluate(A_x)
        assert A0 == A1
        A_x = Polynomial.rand(domain.n)
        A0 = A_x.evaluate_multi(domain.domain)
        A1 = domain.evaluate(A_x)
        assert A0 == A1


def test_omega():
    for i in range(1, 5):
        domain = Domain(domain_config(i))

        A = Polynomial.rand(domain.n)
        A_x = domain.interpolate(A)
        A_wx = domain.omega(A_x)
        Aw = domain.evaluate(A_wx)
        for i in range(A.n()):
            assert A[i] == Aw[i - 1]


def test_mul():
    domain = Domain(domain_config(3))

    n = domain.n
    A = Polynomial.rand(n >> 1)
    B = Polynomial.rand(n >> 1)
    C = Polynomial.rand(n >> 1)
    R0 = A * B
    R1 = domain.mul(A, B)
    assert R0 == R1
    assert domain.div(R1, A) == B
    assert domain.div(R1, B) == A

    A = Polynomial.rand(n)
    B = Polynomial.rand(n)
    C = Polynomial.rand(n)
    assert domain.mul(C, A, B) == domain.mul(C, domain.mul(A, B))
    assert domain.mul(A, B) == domain.mul(B, A)
    assert domain.mul(A, B, C) == domain.mul(C, A, B)

    domain = Domain(domain_config(2))
    A = Polynomial.from_ints([1, 2, 3, 4])
    B = Polynomial.from_ints([2, 4, 6, 8])
    A_x = domain.interpolate(A)
    B_x = domain.interpolate(B)
    C_x = domain.mul(A_x, B_x)
    R0 = domain.evaluate(C_x)
    R1 = A.mul_sample(B)
    assert R0 == R1


def test_domain():
    n = 2
    small_domain = Domain(domain_config(n))
    large_domain = Domain(domain_config(n + 2))

    a = Polynomial.rand(1 << n)
    b = Polynomial.rand(1 << n)
    c = a.mul_sample(b)

    z_x = small_domain.vanishing()
    z_x_dist = z_x.distribute(kx)
    zi = large_domain.evaluate(z_x_dist).inv_sample()

    a_x = small_domain.interpolate(a)
    b_x = small_domain.interpolate(b)
    c_x = small_domain.interpolate(c)
    t0_x = large_domain.mul(a_x, b_x) - c_x
    t0_x = large_domain.coset_div(t0_x, z_x)

    a_x = small_domain.interpolate(a)
    b_x = small_domain.interpolate(b)
    c_x = small_domain.interpolate(c)
    a_x = a_x.distribute(kx)
    b_x = b_x.distribute(kx)
    c_x = c_x.distribute(kx)
    a = large_domain.evaluate(a_x)
    b = large_domain.evaluate(b_x)
    c = large_domain.evaluate(c_x)
    t1 = a.mul_sample(b) - c
    t1 = t1.mul_sample(zi)
    t1_x = large_domain.interpolate(t1)
    t1_x = t1_x.distribute(one / kx)

    assert t1_x == t0_x
    assert t1_x.degree() == a_x.degree() + b_x.degree() - z_x.degree()

    domain = small_domain
    v = [Polynomial.rand(1 << n) for _ in range(10)]
    vi0 = domain.i(*v)
    vi1 = [domain.interpolate(u) for u in v]
    for v0, v1 in zip(vi0, vi1):
        assert v0 == v1


def test_lagrage_polynomial():
    n = 6
    domain = Domain(domain_config(n))
    zeta = Scalar.rand()
    for i in range(domain.n):
        li_x = domain.lagrange_polynomial(i)
        assert li_x(zeta) == domain.lagrange_evaluation(i, zeta)
