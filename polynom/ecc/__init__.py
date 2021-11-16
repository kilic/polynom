from __future__ import annotations

one: Scalar = None
zero: Scalar = None
field_modulus = None
Scalar = None


def init_scalar_field(type):
    global Scalar
    Scalar = type
    global one
    one = Scalar(1)
    global zero
    zero = Scalar(0)


CURVE: Curve = None


class Curve:

    def is_pairing_friendly(self):
        return False

    @staticmethod
    def uncompressed_point_size() -> int:
        pass

    @staticmethod
    def compressed_point_size() -> int:
        pass

    @staticmethod
    def scalar_size() -> int:
        pass

    def from_uncompressed(self, input: bytes):
        pass

    def to_uncompressed(self, p: Point) -> bytes:
        pass

    def scalar_from_bytes(self, input: bytes):
        pass

    def scalar_to_bytes(self, e: Scalar) -> bytes:
        pass


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
    def rand() -> Point:
        return Point.G1(Scalar.rand())

    @staticmethod
    def ZERO():
        return Point(CURVE, CURVE.z1)

    @staticmethod
    def G1(s: Scalar = one):
        p = Point(CURVE, CURVE.g1)
        if s is None:
            return p
        return p * s

    @staticmethod
    def G2(s: Scalar = one):
        p = Point(CURVE, CURVE.g2)
        assert CURVE.is_pairing_friendly()
        if s is None:
            return p
        return p * s

    @staticmethod
    def from_uncompressed(input: bytes) -> Point:
        return Point(CURVE, CURVE.from_uncompressed(input))

    def __init__(self, curve, point):
        self.curve = curve
        self.point = point
        return

    def new(self, point) -> Point:
        return Point(self.curve, point)

    def to_uncompressed(self) -> bytes:
        return self.curve.to_uncompressed(self.point)

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

    def is_on_curve(self):
        return self.curve.is_on_curve_g1(self.point)

    def __eq__(self, other: Point) -> bool:
        return self.curve.eq(self.point, other.point)

    def __repr__(self) -> str:
        return str(self.normalize())
