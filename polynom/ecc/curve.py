from __future__ import annotations
from polynom.scalar import one, Scalar

CURVE: Curve = None


class Curve:

    def is_pairing_friendly(self):
        return False


class PairingFriendlyCurve(Curve):

    def is_pairing_friendly(self) -> bool:
        return True

    def pairing_check(self, pairs: list[tuple[Point, Point]]) -> bool:
        return self._pairing_check(pairs)


def init_ecc(curve: Curve):
    global CURVE
    CURVE = curve


def pairing_check(pairs: list[tuple[Point, Point]]) -> bool:
    assert CURVE.is_pairing_friendly()
    return CURVE.pairing_check(pairs)


class Point():

    @staticmethod
    def ZERO():
        return Point(CURVE, CURVE.z1)

    @staticmethod
    def G1(s: Scalar = one):
        p = Point(CURVE, CURVE.g1)
        return p * s

    @staticmethod
    def G2(s: Scalar = one):
        p = Point(CURVE, CURVE.g2)
        assert CURVE.is_pairing_friendly()
        return p * s

    def __init__(self, curve, point):
        self.curve = curve
        self.point = point
        return

    def new(self, point) -> Point:
        return Point(self.curve, point)

    def __add__(self, other: Point) -> Point:
        point = self.curve.add(self.point, other.point)
        return self.new(point)

    def __sub__(self, other: Point) -> Point:
        neg_other = self.neg(other)
        point = self.curve.add(self.point, neg_other.point)
        return self.new(point)

    def __mul__(self, scalar: Scalar) -> Point:
        new_point = self.curve.mul(self.point, scalar.n)
        return Point(self.curve, new_point)

    def __neg__(self) -> Point:
        new_point = self.curve.neg(self.point)
        return Point(self.curve, new_point)

    def normalize(self):
        return self.curve.normalize(self.point)

    def __eq__(self, other: Point) -> bool:
        return self.curve.add(self.point, other.point)

    def __repr__(self) -> str:
        return str(self.normalize())
