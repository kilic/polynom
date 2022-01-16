from polynom.domain.fft import perform_fft
from typing import Union
from polynom.ecc import Point, Scalar, one, zero
from polynom.utils import log2, pad_scalars, batch_inverse, pad_points
from polynom.polynomial import Polynomial


def calculate_domain(w: Scalar, exp: int, k: int = 1) -> list[Scalar]:
    n = 1 << exp
    W = [Scalar(k)] + (n - 1) * [Scalar(0)]
    for i in range(1, n):
        W[i] = W[i - 1] * w
    return W


class Domain:

    def __init__(self, root_of_unity: Scalar, s: int, exp: int, k: Scalar):

        w = root_of_unity
        w_inv = w.inverse()

        self.exp = exp
        self.s = s
        self.inv_k = one / k
        self.k = k
        self.n = 1 << exp
        self.inv_n = one / self.n

        for _ in range(self.exp, self.s):
            w, w_inv = w**2, w_inv**2

        self.domain = calculate_domain(w, self.exp, one)
        self.inverse_domain = calculate_domain(w_inv, self.exp, one)

        assert self.domain[1] == w
        assert self.inverse_domain[1] == w_inv
        assert batch_inverse(self.domain) == self.inverse_domain

    def extend(self, poly: Polynomial):
        assert poly.n() <= self.n

        k = log2(poly.n())
        u = 1 << self.exp - k
        coeffs: list[Scalar] = []
        for a in poly.coeffs:
            coeffs += [a] + [Scalar(0)] * (u - 1)
        return Polynomial(coeffs)

    def coset(self, k):
        return [k * w for w in self.domain]

    def lagrange_polynomial(self, i: int) -> Polynomial:
        assert i < self.n
        coeffs = [Scalar(0)] * self.n
        coeffs[i] = one
        return self.interpolate(Polynomial(coeffs))

    def lagrange_evaluation_range(self, i: int, j: int, zeta: Scalar) -> Scalar:
        pass

    def lagrange_evaluation(self, i: int, zeta: Scalar) -> Scalar:
        assert i < self.n
        n = self.n
        zeta_n, w = zeta**n, self.domain[i]
        u = (zeta_n - one) * w
        v = (zeta - w) * n
        return u / v

    def new_poly(self, coeffs) -> Polynomial:
        assert len(coeffs) == self.n
        return Polynomial(coeffs)

    def i(self, *input: Union[Polynomial, list[Scalar]]) -> list[Polynomial]:
        return [self.interpolate(poly) for poly in input]

    def interpolate(self, poly: Union[Polynomial, list[Scalar]]) -> Polynomial:
        coeffs = None
        if isinstance(poly, Polynomial):
            assert len(poly) <= self.n
            coeffs = pad_scalars(poly.coeffs, self.n)
            coeffs = [c * self.inv_n for c in perform_fft(coeffs, self.inverse_domain)]
            return Polynomial(coeffs)

        assert isinstance(poly, list)
        coeffs = pad_scalars(poly, self.n)
        coeffs = [c * self.inv_n for c in perform_fft(coeffs, self.inverse_domain)]
        return Polynomial(coeffs)

    def evaluate(self, poly) -> Polynomial:
        coeffs = pad_scalars(poly.coeffs, self.n)
        coeffs = perform_fft(coeffs, self.domain)
        return Polynomial(coeffs)

    def ecc_evaluate(self, points: list[Point]) -> list[Point]:
        assert len(points) <= self.n
        points = pad_points(points, self.n)
        perform_fft(points, self.domain)

    def ecc_interpolate(self, points: list[Point]) -> list[Point]:
        assert len(points) <= self.n
        points = pad_points(points, self.n)
        return [c * self.inv_n for c in perform_fft(points, self.inverse_domain)]

    def w(self) -> Scalar:
        return self.domain[1]

    def w_inv(self) -> Scalar:
        return self.inverse_domain[1]

    def distribute_omega(self, poly: Polynomial) -> Polynomial:
        return poly.distribute(self.w())

    def distribute_zeta(self, poly: Polynomial) -> Polynomial:
        return poly.distribute(self.k)

    def distribute_zeta_inv(self, poly: Polynomial) -> Polynomial:
        return poly.distribute(self.inv_k)

    def vanishing(self) -> Polynomial:
        return Polynomial([-one] + [zero] * (self.n - 1) + [one])

    def mul(self, *v: Polynomial) -> Polynomial:
        assert len(v) > 1
        for u in v:
            assert u.n() <= self.n
            if u.is_zero():
                return Polynomial([zero] * self.n)

        acc = perform_fft(pad_scalars(v[0].coeffs, self.n), self.domain)
        for i in range(1, len(v)):
            v_i_evals = perform_fft(pad_scalars(v[i].coeffs, self.n), self.domain)
            acc = [u * v for u, v in zip(acc, v_i_evals)]

        coeffs = [c * self.inv_n for c in perform_fft(acc, self.inverse_domain)]

        return self.new_poly(coeffs)

    def div(self, a: Polynomial, b: Polynomial) -> Polynomial:
        if a.is_zero() or b.is_zero():
            return Polynomial([zero] * self.n)

        assert a.n() <= self.n
        assert b.n() <= self.n

        a_evals = perform_fft(pad_scalars(a.coeffs, self.n), self.domain)
        b_evals = perform_fft(pad_scalars(b.coeffs, self.n), self.domain)
        b_evals = [1 / u for u in b_evals]

        mul_evals = [u * v for u, v in zip(a_evals, b_evals)]

        coeffs = [c * self.inv_n for c in perform_fft(mul_evals, self.inverse_domain)]

        return self.new_poly(coeffs)

    def coset_div(self, a: Polynomial, b: Polynomial) -> Polynomial:

        if a.is_zero():
            return a.clone()
        u = self.div(self.distribute_zeta(a), self.distribute_zeta(b))
        return self.distribute_zeta_inv(u)
