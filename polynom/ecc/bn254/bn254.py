from py_ecc.optimized_bn128 import pairing, add, multiply, G1, G2, Z1, Z2, FQ12, normalize, neg, eq, is_on_curve, b, b2, FQ
from polynom.ecc.bn254.scalar import Scalar
from polynom.ecc import PairingFriendlyCurve, Point


class bn254(PairingFriendlyCurve):

    @staticmethod
    def uncompressed_point_size() -> int:
        return 64

    @staticmethod
    def compressed_point_size() -> int:
        return 32

    @staticmethod
    def scalar_size() -> int:
        return 32

    @staticmethod
    def rom_scalar_size() -> int:
        return 64

    def from_uncompressed(self, input: bytes):
        u = self.uncompressed_point_size()
        assert len(input) >= u
        u = u // 2
        x = int.from_bytes(input[:u], "little")
        y = int.from_bytes(input[u:u * 2], "little")
        point = (FQ(x), FQ(y), FQ(1))
        assert self.is_on_curve_g1(point)
        return point

    def to_uncompressed(self, p: Point) -> bytes:
        normalized = normalize(p)
        x = normalized[0]
        y = normalized[1]
        return self.scalar_to_bytes(x) + self.scalar_to_bytes(y)

    def scalar_from_bytes(self, input: bytes):
        return Scalar.from_32(input)

    def scalar_to_bytes(self, e) -> bytes:
        return e.n.to_bytes(bn254.scalar_size(), 'little')

    def __init__(self):
        self.add = add
        self.mul = multiply
        self.neg = neg
        self.eq = eq
        self.normalize = normalize
        self.z1 = Z1
        self.g1 = G1

        self.g2 = G2
        self.z2 = Z2

    def _pairing_check(self, pairs) -> bool:
        acc = FQ12.one()
        # TODO: use multi miller look
        for pair in pairs:
            acc = acc * pairing(pair[0].point, pair[1].point)

        return acc == FQ12.one()

    def is_on_curve_g1(self, p: Point) -> bool:
        return is_on_curve(p, b)


BN254 = bn254()