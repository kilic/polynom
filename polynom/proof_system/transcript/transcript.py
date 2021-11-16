from polynom.ecc import Point, Scalar, CURVE
from polynom.proof_system.transcript.hasher import Hasher


class Transcript():

    def __init__(self, hasher):
        self.hasher = hasher

    def challenge(self,) -> Scalar:
        return self.hasher.challenge()

    def write_point_to_state(self, point: Point) -> bytes:
        return self.hasher.update_point(point)

    def write_scalar_to_state(self, e: Scalar) -> bytes:
        return self.hasher.update_scalar(e)


class TranscriptRead(Transcript):

    def __init__(self, hasher: Hasher, message: bytes):
        hasher.clean_state()
        super().__init__(hasher)
        self.message = message
        self.offset = 0

    def read_4_bytes(self) -> int:
        ret = int.from_bytes(self.message[self.offset:self.offset + 4], "little")
        self.offset += 4
        return ret

    def read_point(self) -> Point:
        point_size = CURVE.uncompressed_point_size()
        assert len(self.message) >= self.offset + point_size
        point = Point.from_uncompressed(self.message[self.offset:self.offset + point_size])
        self.offset += point_size
        self.write_point_to_state(point)
        return point

    def read_scalar(self) -> Scalar:
        scalar_size = CURVE.scalar_size()
        assert len(self.message) >= self.offset + scalar_size
        scalar = CURVE.scalar_from_bytes(self.message[self.offset:self.offset + scalar_size])
        self.offset += scalar_size
        self.write_scalar_to_state(scalar)
        return scalar


class TranscriptWrite(Transcript):

    def __init__(self, hasher: Hasher):
        hasher.clean_state()
        super().__init__(hasher)
        self.message = bytes()

    def write_4_bytes(self, n: int):
        self.message += n.to_bytes(4, "little")

    def write_point(self, point: Point):
        self.message += self.write_point_to_state(point)

    def write_scalar(self, e: Scalar):
        self.message += self.write_scalar_to_state(e)

    def get_message(self) -> bytes:
        return self.message
