from py_ecc.optimized_bn128 import pairing, add, multiply, G1, G2, Z1, Z2, FQ12, normalize, neg, eq
from polynom.ecc.curve import PairingFriendlyCurve


class bn254(PairingFriendlyCurve):

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


BN254 = bn254()