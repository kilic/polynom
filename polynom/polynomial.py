from __future__ import annotations
from typing import Tuple
from polynom.ecc import Scalar, zero, one
from polynom.utils import pad_scalars, trim_zeros


def evaluate(point: Scalar, *inputs: Polynomial) -> list[Scalar]:
    return [input(point) for input in inputs]


def lagrange_interpolation(xy: list[Tuple(Scalar, Scalar)]) -> Polynomial:

    L_x = Polynomial.zero()
    size = len(xy)
    for i in range(size):
        l_x = Polynomial.one()
        denom = Scalar(1)
        x_i, y_i = xy[i]
        for j in range(size):
            x_j, _ = xy[j]
            if j == i:
                continue
            contrib = Polynomial.degree_one(x_j)
            l_x = l_x * contrib
            denom = denom * (x_i - x_j)
        l_x = l_x * (y_i / denom)
        L_x = L_x + l_x

    return L_x


def barycentric_preprocess(points: list[Scalar]) -> list[tuple[Scalar, Scalar]]:

    weights = []
    size = len(points)
    for i in range(size):
        acc = one
        x_i = points[i]
        for j in range(size):
            if i == j:
                continue
            x_j = points[j]
            dif = (x_i - x_j)
            # expect distinct points
            assert dif != zero
            acc = acc * dif
        # one can apply batch inversion to make it more effective
        weights.append((x_i, one / acc))
    return weights


def barycentric_evaluation(weights: list[tuple[Scalar, Scalar]], evaluations: list[Scalar], z: Scalar):
    size = len(evaluations)
    assert len(weights) == size
    # one can apply batch inversion to make it more effective
    denom_coeffs = [weights[i][1] * (one / (z - weights[i][0])) for i in range(size)]
    num_coeffs = [coeff * e for coeff, e in zip(denom_coeffs, evaluations)]
    return sum(num_coeffs) / sum(denom_coeffs)


class Polynomial:

    @staticmethod
    def from_ints(v: list[int]):

        return Polynomial([Scalar(e) for e in v])

    @staticmethod
    def degree_one(z: Scalar):
        return Polynomial([-z, Scalar(1)])

    @staticmethod
    def zero():

        return Polynomial([])

    @staticmethod
    def one():

        return Polynomial([Scalar(1)])

    @staticmethod
    def rand(n: int):

        ret = Polynomial([Scalar.rand() for _ in range(n)])
        assert not ret.is_zero()
        return ret

    def __init__(self, coeffs: list[Scalar]):

        self.coeffs = coeffs
        self.evals = {}

    def __call__(self, z) -> Scalar:
        eval = self.evaluate(z)
        assert isinstance(eval, Scalar)
        return eval

    def debug_str(self, label, verbose=False) -> str:
        res = label
        res += "\tlen: " + str(self.n())
        res += "\tdegree: " + str(self.degree())
        if verbose:
            res = res + "\n"
            for c in self.coeffs:
                res = res + hex(c.n) + "\n"
        return res

    def debug(self, label="") -> Polynomial:
        print(self.debug_str(label, False))
        return self

    def debug_verbose(self, label="") -> Polynomial:
        print(self.debug_str(label, True))
        return self

    def clone(self) -> Polynomial:
        return Polynomial(self.coeffs[:])

    def pad(self, n):
        coeffs = pad_scalars(self.coeffs, n)
        return Polynomial(coeffs)

    def eq(self, other) -> bool:
        a = self.trim_zeros()
        b = other.trim_zeros()
        if a.n() != b.n():
            return False
        for ai, bi in zip(a.coeffs, b.coeffs):
            if ai != bi:
                return False
        return True

    def n(self) -> int:
        return len(self.coeffs)

    def is_zero(self):
        for e in self.coeffs:
            if e != zero:
                return False
        return True

    def degree(self) -> int:
        res = len(self.coeffs) - 1
        for c in reversed(self.coeffs):
            if c == zero:
                res -= 1
            else:
                return res
        return res

    def trim_zeros(self) -> Polynomial:
        return Polynomial(trim_zeros(self.coeffs))

    def evaluate(self, x: Scalar) -> Scalar:
        if x.n in self.evals:
            return self.evals[x.n]
        acc = x.zero()
        for a in reversed(self.coeffs):
            acc = x * acc + a
        self.evals[x.n] = acc
        return acc

    def evaluate_multi(self, xs: list[Scalar]) -> Polynomial:
        samples: list[Scalar] = []
        for x in xs:
            samples.append(self.evaluate(x))
        return Polynomial(samples)

    def add(self, other) -> Polynomial:

        if isinstance(other, Scalar):
            coeffs = self.coeffs[:]
            coeffs[0] = coeffs[0] + other
            return Polynomial(coeffs)

        if isinstance(other, list):
            other = Polynomial(other)

        if self.is_zero():
            return other.clone()
        if other.is_zero():
            return self.clone()

        n = max(other.n(), self.n())
        coeffs = [self[i] + other[i] for i in range(n)]
        return Polynomial(coeffs)

    def sub(self, other) -> Polynomial:

        if isinstance(other, Scalar):
            coeffs = self.coeffs[:]
            coeffs[0] = coeffs[0] - other
            return Polynomial(coeffs)

        if isinstance(other, list):
            other = Polynomial(other)

        if self.is_zero():
            return other.clone()
        if other.is_zero():
            return self.clone()

        n = max(other.n(), self.n())
        coeffs = [self[i] - other[i] for i in range(n)]
        return Polynomial(coeffs)

    def neg(self) -> Polynomial:
        return Polynomial([-a for a in self.coeffs])

    def scale(self, k) -> Polynomial:
        return Polynomial([a * k for a in self.coeffs])

    def distribute(self, k) -> Polynomial:
        coeffs = []
        acc = 1
        for a in self.coeffs:
            coeffs.append(a * acc)
            acc = acc * k
        return Polynomial(coeffs)

    def mul_naive(self, b) -> Polynomial:
        a = self
        if isinstance(b, list):
            b = Polynomial(b)
        c = (a.n() + b.n() - 1) * [Scalar(0)]
        for i, u in enumerate(a.coeffs):
            for j, v in enumerate(b.coeffs):
                c[i + j] += u * v
        return Polynomial(c)

    def mul_sample(self, other) -> Polynomial:
        n = min(self.n(), other.n())
        return Polynomial([self.coeffs[i] * other.coeffs[i] for i in range(n)])

    def inv_sample(self) -> Polynomial:
        # TODO: use batch inversion
        return Polynomial([one / a for a in self.coeffs])

    def __eq__(self, other):
        return self.eq(other)

    def __ne__(self, other):
        return not self.eq(other)

    def __add__(self, other):
        return self.add(other)

    def __neg__(self):
        return self.neg()

    def __sub__(self, other):
        return self.sub(other)

    def __mul__(self, other):
        if isinstance(other, Scalar):
            return self.scale(other)
        return self.mul_naive(other)

    def __repr__(self) -> str:
        return self.debug_str("", False)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return Polynomial(self.coeffs[i][:])
        return self.coeffs[i] if self.n() > i else zero

    def __len__(self):
        return len(self.coeffs)

    # def div_vanishing(self, z) -> Polynomial:
    #     if self.is_zero():
    #         return self.clone()
    #     u = self.distribute(k) / z.distribute(k)
    #     return u.distribute(one / k)
