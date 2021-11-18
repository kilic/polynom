from __future__ import annotations
from polynom.commitment.bdfg.prover import BDFGProver
from polynom.commitment.bdfg.verifier import BDFGVerifier
from polynom.commitment.gwc.gwc import GWCProver, GWCVerifier
from polynom.commitment.kzg import KZGProver, KZGVerifier
from polynom.ecc import Point
from polynom.ecc import Scalar
from polynom.domain import Domain
from polynom.proof_system.transcript.hasher import Hasher

ZERO_POINT = Point.ZERO()


class KZGSetup:

    def __init__(self, bases: list[Point], domain: Domain, X_2, hasher: Hasher):
        self.bases = bases
        self.domain = domain
        self.X_2 = X_2
        self.G_1 = Point.G1()
        self.nG2 = -Point.G2()
        self.hasher = hasher

    @staticmethod
    def new(domain: Domain, hasher: Hasher) -> KZGSetup:
        # FIX: only for bn254???
        s = Scalar(0x0330fa29c0b79377aa26b9f89ad6c94912201b4c8a854c4fe1db0aae5d3e3139)
        n = domain.n
        bases = [Point.G1()]
        r = Scalar(1)
        for _ in range(n - 1):
            r = r * s
            bases.append(Point.G1(r))
        sG = Point.G2(s)
        return KZGSetup(bases, domain, sG, hasher)

    def prover_kzg(self) -> KZGProver:
        return KZGProver(self.hasher, self.bases, self.domain)

    def verifier_kzg(self) -> KZGVerifier:
        return KZGVerifier(self.hasher, self.G_1, self.X_2, self.domain.w())

    def prover_gwc(self) -> GWCProver:
        return GWCProver(self.hasher, self.bases, self.domain)

    def verifier_gwc(self) -> GWCVerifier:
        return GWCVerifier(self.hasher, self.G_1, self.X_2, self.domain.w())

    def prover_bdfg(self) -> BDFGProver:
        return BDFGProver(self.hasher, self.bases, self.domain)

    def verifier_bdfg(self) -> BDFGVerifier:
        return BDFGVerifier(self.hasher, self.G_1, self.X_2, self.domain.w())
