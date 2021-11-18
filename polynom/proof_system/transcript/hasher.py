import hashlib
from polynom.ecc import Point, CURVE, Scalar


class Hasher:

    def __init__(self, person: bytes, scalar_prefix: bytes, challenge_prefix: bytes, point_prefix: bytes):
        pass

    def update_scalar(self, e: Scalar):
        pass

    def update_point(self, e: Point):
        pass

    def challenge(self) -> Scalar:
        pass

    def clean_state(self):
        pass


from hashlib import sha256


class SHA256(Hasher):

    @staticmethod
    def new_state(person: bytes):
        hasher = sha256()
        hasher.update(person)
        return hasher

    def __init__(self, person: bytes, scalar_prefix: bytes, challenge_prefix: bytes, point_prefix: bytes):
        self.scalar_prefix = scalar_prefix
        self.point_prefix = point_prefix
        self.challenge_prefix = challenge_prefix
        self.person = person
        self.hasher = self.clean_state()

    def update_scalar(self, e: Scalar) -> bytes:
        self.hasher.update(self.scalar_prefix)
        in_bytes = CURVE.scalar_to_bytes(e)
        self.hasher.update(in_bytes)
        return in_bytes

    def update_point(self, e: Point) -> bytes:
        self.hasher.update(self.point_prefix)
        in_bytes = e.to_uncompressed()
        self.hasher.update(in_bytes)
        return in_bytes

    def challenge(self) -> Scalar:
        # FIX: this is not secure way to squeeze challenge
        # We must derive a scalar from larger input
        self.hasher.update(self.challenge_prefix + (0).to_bytes(1, 'little'))
        u_0 = self.hasher.digest()
        return Scalar.from_32(u_0)
        # self.hasher.update(self.challenge_prefix + (1).to_bytes(1, 'little'))
        # u_1 = self.hasher.digest()
        # return Scalar.from_64(u_0 + u_1)

    def clean_state(self):
        self.hasher = SHA256.new_state(self.person)


class Keccak256(Hasher):
    pass


# >>> import hashlib
# >>> m = hashlib.sha256()
# >>> m.update(b"Nobody inspects")
# >>> m.update(b" the spammish repetition")
# >>> m.digest()
# b'\x03\x1e\xdd}Ae\x15\x93\xc5\xfe\\\x00o\xa5u+7\xfd\xdf\xf7\xbcN\x84:\xa6\xaf\x0c\x95\x0fK\x94\x06'
# >>> m.digest_size
# 32
# >>> m.block_size
