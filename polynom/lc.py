from polynom.scalar import Scalar, zero, one
from polynom.ecc.curve import Point
from polynom.polynomial import Polynomial

ZERO_POINT = Point.ZERO()


class LinearCombination:

    def __init__(self, e: Scalar):
        self.e = e

    def combine_poly(self, *coeffs: Polynomial) -> Polynomial:
        acc, e = Polynomial([zero]), one
        for c in coeffs:
            acc = acc + c * e
            e = e * self.e
        return acc

    def combine_fr(self, *coeffs: Scalar) -> Scalar:
        acc, e = zero, one
        for c in coeffs:
            acc = acc + c * e
            e = e * self.e
        return acc

    def combine_points_for_degree(self, degree: int, *input: tuple[Point, Scalar]) -> Point:
        assert len(input) > 0
        e, acc = self.e**degree, ZERO_POINT
        for (point, scalar) in input:
            acc = acc + point * (scalar * e)
        return acc

    def multiexp_with_aux(self, degree: int, *input: tuple[Point, Scalar]) -> Point:
        assert len(input) > 0
        e, acc = self.e**degree, ZERO_POINT
        for (point, scalar) in input:
            acc = acc + point * (scalar * e)
        return acc

    def combine_points(self, *points: Point) -> Point:
        acc, e = ZERO_POINT, one
        for point in points:
            acc = acc + point * e
            e = e * self.e
        return acc

    def combine_ecc_with_aux(self, *inputs: tuple[Point, Scalar]) -> tuple[Point, Point]:
        accW, accR, e = ZERO_POINT, ZERO_POINT, one
        for (point, zeta) in inputs:
            accW = accW + point * e
            accR = accR + point * (e * zeta)
            e = e * self.e
        return accW, accR
