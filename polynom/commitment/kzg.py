from __future__ import annotations
from polynom.ecc.curve import Point, pairing_check
from polynom.scalar import Scalar, one, zero
from polynom.polynomial import Polynomial
from polynom.domain import Domain
from polynom.lc import LinearCombination

ZERO_POINT = Point.ZERO()


def kzg_setup(n, domain: Domain) -> tuple[KZGProver, KZGVerifier]:
    s = Scalar(0x0330fa29c0b79377aa26b9f89ad6c94912201b4c8a854c4fe1db0aae5d3e3139)
    bases = [Point.G1()]
    r = Scalar(1)
    for _ in range(n - 1):
        r = r * s
        bases.append(Point.G1(r))
    sG = Point.G2(s)
    return KZGProver(bases, domain), KZGVerifier(sG)


class KZGProof:

    z: Scalar
    W: Point
    F: Point
    e: Scalar

    def __init__(
        self,
        z: Scalar,
        W: Point,
        F: Point,
        e: Scalar,
    ) -> None:
        self.z = z
        self.W = W
        self.F = F
        self.e = e


class KZGProver:

    def __init__(self, bases: list[Point], domain: Domain):
        self.bases = bases
        self.domain = domain

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

    def witness(self, z: Scalar, p_x: Polynomial) -> Point:
        y = p_x(z)
        Q_x = p_x - y
        root_z = Polynomial([-z, Scalar(1)])
        W_x = self.domain.div(Q_x, root_z)
        return self.commit(W_x)

    def witness_single_batch(self, challenge: LinearCombination, z: Scalar, _polys: list[Polynomial]) -> Point:
        polys = [p_x - p_x(z) for p_x in _polys]
        u_x = challenge.combine_poly(*polys)
        root_z = Polynomial([-z, one])
        w_x = self.domain.div(u_x, root_z)
        return self.commit(w_x)

    def witness_batch(self, inputs: list[tuple[LinearCombination, Scalar, list[Polynomial]]]) -> list[Point]:
        W: list[Point] = []
        for input in inputs:
            challange, z, polys = input
            W.append(self.witness_single_batch(challange, z, polys))
        return W


class KZGVerifier:

    def __init__(self, sG2: Point):
        self.sG2 = sG2
        self.nG2 = -Point.G2()

    def verify(
        self,
        proof: KZGProof,
    ):
        E = Point.G1(-proof.e)
        zW = proof.W * proof.z

        pairs = [(self.sG2, proof.W), (self.nG2, zW + proof.F + E)]
        return pairing_check(pairs)

    def verify_batch(self, multipoint_eval_challange: LinearCombination, proofs: list[KZGProof]):
        z_W_i: list[tuple[Point, Scalar]] = []
        e: list[Scalar] = []
        F: list[Point] = []
        for proof in proofs:
            z = proof.z
            F.append(proof.F)
            e.append(proof.e)

            z_W_i.append((proof.W, z))

        F_combined = multipoint_eval_challange.combine_points(*F)
        e_combined = multipoint_eval_challange.combine_fr(*e)
        E = Point.G1(-e_combined)
        W, zW = multipoint_eval_challange.combine_ecc_with_aux(*z_W_i)

        pairs = [(self.sG2, W), (self.nG2, zW + E + F_combined)]
        return pairing_check(pairs)
