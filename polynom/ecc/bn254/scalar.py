from py_ecc.fields.field_elements import FQ
import random

MODULUS = (21888242871839275222246405745257275088548364400416034343698204186575808495617)
T = 1 << 256
R = T % MODULUS
R2 = (R * R) % MODULUS
R3 = (R2 * R) % MODULUS


class Scalar(FQ):

    field_modulus = MODULUS

    @staticmethod
    def rand(n: int = 1):
        if n == 1:
            return Scalar(random.randint(0, MODULUS))
        return [Scalar.rand() for i in range(n)]

    @staticmethod
    def one():
        return Scalar(1)

    @staticmethod
    def zero():
        return Scalar(0)

    @staticmethod
    def from_64(input):
        assert len(input) >= 64
        # return FQ(int.from_bytes(input[:u], "little"))
        assert False  # TODO: implement

    @staticmethod
    def from_32(input):
        assert len(input) >= 32
        return Scalar(int.from_bytes(input[:32], "little"))

    def to_bytez(self) -> bytes:
        return self.n.to_bytes(32, "little")

    def __repr__(self) -> str:
        return hex(self.n)

    def __hash__(self) -> int:
        return hash(self.n)

    def __lt__(self, other) -> bool:
        return self.n < other.n

    def __le__(self, other) -> bool:
        return self.n <= other.n

    def __gt__(self, other) -> bool:
        return self.n > other.n

    def __ge__(self, other) -> bool:
        return self.n >= other.n
