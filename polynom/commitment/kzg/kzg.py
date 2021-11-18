from __future__ import annotations
from polynom.commitment.kzg_base import KZGProverBase, KZGVerifierBase
from polynom.ecc import Point, pairing_check
from polynom.ecc import Scalar, one
from polynom.polynomial import Polynomial, evaluate
from polynom.lc import LinearCombination


class KZGProver(KZGProverBase):

    def create_proof(self, P_x: Polynomial) -> bytes:

        transcript = self.new_transcript()

        P = self.commit(P_x)
        transcript.write_point(P)
        z = transcript.challenge()

        eval = P_x(z)
        Q_x = P_x - eval
        root_z = Polynomial([-z, Scalar(1)])
        W_x = self.domain.div(Q_x, root_z)
        transcript.write_scalar(eval)

        W = self.commit(W_x)
        transcript.write_point(W)
        return transcript.get_message()

    def create_proof_batch(self, polys: list[Polynomial]) -> Point:

        n = len(polys)
        assert n > 1

        transcript = self.new_transcript()

        commitments = self.c(*polys)
        for c in commitments:
            transcript.write_point(c)

        z = transcript.challenge()
        alpha = LinearCombination(transcript.challenge())

        evals = evaluate(z, *polys)
        for e in evals:
            transcript.write_scalar(e)

        polys = [poly - eval for eval, poly in zip(evals, polys)]
        u_x = alpha.combine_poly(*polys)
        root_z = Polynomial([-z, one])
        w_x = self.domain.div(u_x, root_z)

        W = self.commit(w_x)
        transcript.write_point(W)
        return transcript.get_message()


class KZGVerifier(KZGVerifierBase):

    def verify(self, proof: bytes) -> bool:

        transcript = self.new_transcript(proof)

        F = transcript.read_point()
        z = transcript.challenge()
        eval = transcript.read_scalar()
        W = transcript.read_point()

        E = Point.G1(-eval)

        zW = W * z
        pairs = [(self.X_2, W), (self.n_G_2, zW + F + E)]

        return pairing_check(pairs)

    def verify_batch(self, commitment_size: int, proof: bytes) -> bool:

        transcript = self.new_transcript(proof)

        commitments = [transcript.read_point() for _ in range(commitment_size)]

        z = transcript.challenge()
        alpha = LinearCombination(transcript.challenge())

        F = alpha.combine_points(*commitments)
        evals = [transcript.read_scalar() for _ in range(commitment_size)]
        e_combined = alpha.combine_fr(*evals)
        E = Point.G1(-e_combined)

        W = transcript.read_point()
        z_W = W * z

        pairs = [(self.X_2, W), (self.n_G_2, z_W + F + E)]

        return pairing_check(pairs)
