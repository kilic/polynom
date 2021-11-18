from __future__ import annotations
from polynom.commitment.kzg_base import KZGProverBase, KZGVerifierBase
from polynom.ecc import Point, pairing_check
from polynom.ecc import one
from polynom.polynomial import Polynomial, evaluate
from polynom.lc import LinearCombination


class GWCKey():

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


class GWCProver(KZGProverBase):

    def create_proof(self, polys: list[Polynomial], key: GWCKey) -> bytes:

        commitment_size = key.commitment_size()
        assert commitment_size > 0
        assert commitment_size == len(polys)

        transcript = self.new_transcript()

        [transcript.write_point(c) for c in self.c(*polys)]

        z = transcript.challenge()

        for shift_val in key.shift_values():

            eval_point = z * (self.domain.w()**shift_val)
            polys_to_eval = [polys[i] for i in key.poly_indexes(shift_val)]
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


class GWCVerifier(KZGVerifierBase):

    def verify(self, key: GWCKey, proof: bytes) -> bool:

        commitment_size = key.commitment_size()
        assert commitment_size > 0

        transcript = self.new_transcript(proof)

        commitments = [transcript.read_point() for _ in range(commitment_size)]

        z = transcript.challenge()
        witnesses, witnesses_mul_evals, combined_commitments, combined_evals = [], [], [], []

        for shift_val in key.shift_values():

            eval_point = z * (self.w**shift_val)
            commitments_to_open = [commitments[i] for i in key.poly_indexes(shift_val)]
            evals = [transcript.read_scalar() for _ in range(len(commitments_to_open))]

            alpha = LinearCombination(transcript.challenge())

            W = transcript.read_point()

            combined_commitments.append(alpha.combine_points(*commitments_to_open))
            combined_evals.append(alpha.combine_fr(*evals))
            witnesses.append(W)
            witnesses_mul_evals.append(W * eval_point)

        multi_open_challenge = LinearCombination(transcript.challenge())

        W = multi_open_challenge.combine_points(*witnesses)
        z_W = multi_open_challenge.combine_points(*witnesses_mul_evals)
        combined_evals = multi_open_challenge.combine_fr(*combined_evals)
        E = Point.G1(-combined_evals)
        F = multi_open_challenge.combine_points(*combined_commitments)

        pairs = [(self.X_2, W), (self.n_G_2, z_W + F + E)]

        return pairing_check(pairs)