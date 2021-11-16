from __future__ import annotations
from polynom.ecc import Point, pairing_check
from polynom.ecc import Scalar, one, zero
from polynom.polynomial import Polynomial, evaluate
from polynom.domain import Domain
from polynom.lc import LinearCombination
from polynom.proof_system.transcript.hasher import Hasher
from polynom.proof_system.transcript.transcript import TranscriptRead, TranscriptWrite

ZERO_POINT = Point.ZERO()


class KZGSetup:

    def __init__(self, bases: list[Point], domain: Domain, sG2, hasher: Hasher):
        self.bases = bases
        self.domain = domain
        self.sG2 = sG2
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

    def prover(self) -> KZGProver:
        return KZGProver(self.hasher, self.bases, self.domain)

    def verifier(self) -> KZGVerifier:
        return KZGVerifier(self.hasher, self.sG2, self.domain.w())


class BatchGWC():

    @staticmethod
    def empty() -> BatchGWC:
        return BatchGWC({})

    def __init__(self, map: dict[int, list[int]]):
        self.map = map

    def add_to_map(self, shift_val: int, poly_index: int):
        if shift_val not in self.map:
            self.map[shift_val] = []
        self.map[shift_val].append(poly_index)

    def shift_values(self) -> list[int]:
        return self.map.keys()

    def poly_indexes(self, shift_val) -> list[int]:
        return self.map[shift_val]

    def set_size(self) -> int:
        return len(self.map.keys)

    def commitment_size(self) -> int:
        max_index = 0
        indexes = {}
        for value in self.map.values():
            for index in value:
                max_index = max(index, max_index)
                indexes[index] = True
        size = len(indexes.keys())
        assert size - 1 == max_index
        return size


class KZGProver():

    def __init__(self, hasher: Hasher, bases: list[Point], domain: Domain):
        self.bases = bases
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

    def new_transcript(self) -> TranscriptWrite:
        return TranscriptWrite(self.hasher)

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

    def create_proof_single_point_open_batch(self, polys: list[Polynomial]) -> Point:

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

    def create_proof_multi_point_open_batch_gwc(self, polys: list[Polynomial], context: BatchGWC) -> bytes:

        commitment_size = context.commitment_size()
        assert commitment_size > 0
        assert commitment_size == len(polys)

        transcript = self.new_transcript()

        [transcript.write_point(c) for c in self.c(*polys)]

        z = transcript.challenge()

        for shift_val in context.shift_values():

            eval_point = z * (self.domain.w()**shift_val)
            polys_to_eval = [polys[i] for i in context.poly_indexes(shift_val)]
            evals = evaluate(eval_point, *polys_to_eval)

            [transcript.write_scalar(e) for e in evals]

            alpha = LinearCombination(transcript.challenge())

            polys_w_root = [poly - eval for eval, poly in zip(evals, polys_to_eval)]
            u_x = alpha.combine_poly(*polys_w_root)
            root_z = Polynomial([-eval_point, one])
            w_x = self.domain.div(u_x, root_z)

            W = self.commit(w_x)
            transcript.write_point(W)

        return transcript.get_message()


class KZGVerifier():

    def __init__(self, hasher: Hasher, sG2: Point, w: Scalar):
        self.hasher = hasher
        self.sG2 = sG2
        self.nG2 = -Point.G2()
        self.w = w

    def new_transcript(self, proof: bytes) -> TranscriptRead:
        return TranscriptRead(self.hasher, proof)

    def verify(self, proof: bytes) -> bool:

        transcript = self.new_transcript(proof)

        F = transcript.read_point()
        z = transcript.challenge()
        eval = transcript.read_scalar()
        W = transcript.read_point()

        E = Point.G1(-eval)

        zW = W * z
        pairs = [(self.sG2, W), (self.nG2, zW + F + E)]

        return pairing_check(pairs)

    def verify_single_point_open_batch(self, commitment_size: int, proof: bytes) -> bool:

        transcript = self.new_transcript(proof)

        commitments = [transcript.read_point() for _ in range(commitment_size)]

        z = transcript.challenge()
        alpha = LinearCombination(transcript.challenge())

        F = alpha.combine_points(*commitments)
        evals = [transcript.read_scalar() for _ in range(commitment_size)]
        e_combined = alpha.combine_fr(*evals)
        E = Point.G1(-e_combined)

        W = transcript.read_point()
        zW = W * z

        pairs = [(self.sG2, W), (self.nG2, zW + F + E)]

        return pairing_check(pairs)

    def verify_multi_point_open_batch_gwc(self, context: BatchGWC, proof: bytes) -> bool:

        commitment_size = context.commitment_size()
        assert commitment_size > 0

        transcript = self.new_transcript(proof)

        commitments = [transcript.read_point() for _ in range(commitment_size)]

        z = transcript.challenge()
        witnesses, witnesses_mul_evals, combined_commitments, combined_evals = [], [], [], []

        for shift_val in context.shift_values():

            eval_point = z * (self.w**shift_val)
            commitments_to_open = [commitments[i] for i in context.poly_indexes(shift_val)]
            evals = [transcript.read_scalar() for _ in range(len(commitments_to_open))]

            alpha = LinearCombination(transcript.challenge())

            W = transcript.read_point()

            combined_commitments.append(alpha.combine_points(*commitments_to_open))
            combined_evals.append(alpha.combine_fr(*evals))
            witnesses.append(W)
            witnesses_mul_evals.append(W * eval_point)

        multi_open_challenge = LinearCombination(transcript.challenge())

        W = multi_open_challenge.combine_points(*witnesses)
        zW = multi_open_challenge.combine_points(*witnesses_mul_evals)
        combined_evals = multi_open_challenge.combine_fr(*combined_evals)
        E = Point.G1(-combined_evals)
        F = multi_open_challenge.combine_points(*combined_commitments)

        pairs = [(self.sG2, W), (self.nG2, zW + F + E)]

        return pairing_check(pairs)
