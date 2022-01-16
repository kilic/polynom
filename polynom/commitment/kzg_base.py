from __future__ import annotations
from polynom.ecc import Point, Scalar
from polynom.polynomial import Polynomial
from polynom.domain import Domain
from polynom.proof_system.transcript.hasher import Hasher
from polynom.proof_system.transcript.transcript import TranscriptRead, TranscriptWrite

ZERO_POINT = Point.ZERO()


class KZGProverBase():

    def __init__(self, hasher: Hasher, bases: list[Point], inverse_bases: list[Point], domain: Domain):

        self.bases = bases
        self.inverse_bases = inverse_bases
        self.domain = domain
        self.hasher = hasher

    def n(self):

        return len(self.bases)

    def c(self, *inputs: Polynomial) -> list[Point]:

        return [self.commit(p_x) for p_x in inputs]

    def commit(self, p_x: Polynomial) -> Point:

        assert self.n() >= p_x.n()
        acc = ZERO_POINT
        for i, a in enumerate(p_x.coeffs):
            acc = acc + self.bases[i] * a
        return acc

    def commit_lagrange(self, p_x: Polynomial) -> Point:

        assert self.n() >= p_x.n()
        acc = ZERO_POINT
        for i, a in enumerate(p_x.coeffs):
            acc = acc + self.inverse_bases[i] * a
        return acc

    def new_transcript(self) -> TranscriptWrite:

        return TranscriptWrite(self.hasher)


class KZGVerifierBase():

    def __init__(self, hasher: Hasher, G: Point, X_2: Point, w: Scalar):

        self.hasher = hasher
        self.G = G
        self.X_2 = X_2
        self.n_G_2 = -Point.G2()
        self.w = w

    def new_transcript(self, proof: bytes) -> TranscriptRead:

        return TranscriptRead(self.hasher, proof)
