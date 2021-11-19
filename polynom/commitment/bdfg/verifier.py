from __future__ import annotations, barry_as_FLUFL
from polynom.commitment.kzg_base import KZGVerifierBase
from polynom.commitment.bdfg.common import BatchBDFGCommon, MultiBDFGCommon, vanising_at
from polynom.ecc import Scalar, pairing_check
from polynom.lc import LinearCombination
from polynom.polynomial import Polynomial, lagrange_interpolation


class MultiBDFGVerifierKey(MultiBDFGCommon):

    def low_degree_equivalent(self, evals: list[Scalar], z: Scalar):

        eval_points = self.eval_points(z)
        assert len(eval_points) == len(evals)
        points = [(x, y) for x, y in zip(eval_points, evals)]
        return lagrange_interpolation(points)


class BatchBDFGVerifierKey(BatchBDFGCommon):

    def __init__(self, w: Scalar, openings: list[MultiBDFGVerifierKey]):
        self.openings = openings
        self.w = w

    def inverse_vanishing(self, multi_open_index: int, z: Scalar) -> Polynomial:
        # TODO: double check the order
        inv_eval_points = self.eval_points(z) ^ set(self.openings[multi_open_index].eval_points(z))
        return vanising_at(inv_eval_points)

    def opening_size(self) -> int:
        return len(self.openings)


class BDFGVerifier(KZGVerifierBase):

    def new_multi_key(self, shifts: list[int]) -> MultiBDFGVerifierKey:
        return MultiBDFGVerifierKey(self.w, shifts)

    def new_batch_key(self, openings: list[MultiBDFGVerifierKey]) -> BatchBDFGVerifierKey:
        return BatchBDFGVerifierKey(self.w, openings)

    def verify_single(self, proof: bytes, key: MultiBDFGVerifierKey) -> bool:

        transcript = self.new_transcript(proof)

        # read the comitment
        F = transcript.read_point()
        # get evaluation seed
        z = transcript.challenge()

        # read evaluations
        evals = [transcript.read_scalar() for _ in range(key.opening_size())]

        # read first witness points
        W = transcript.read_point()

        # get linearisation point
        x = transcript.challenge()

        # read second witness points
        W_2 = transcript.read_point()

        # construct r(X) that is low degree equivalent of `f(X)`
        # where `f(set(z)) = r(set(z))`
        # then make a soft commitment to it `R = r(z) * G`
        r_x = key.low_degree_equivalent(evals, z)
        R = self.G * r_x(x)

        # pairing check:

        # `e(W', X) == e(z_x(z) * W' + L, G)``

        # first we reconstruct linearisation point `L`
        # `L = com(f(x)) - R - com(h(x)) * Z(x)``
        z_x = key.vanising(z)
        L = F - R - (W * z_x(x))

        # calcualate the escape from G2 term `z_x(z) * W'
        x_W_2 = W_2 * x

        # finally apply the pairing
        pairs = [(self.X_2, W_2), (self.n_G_2, x_W_2 + L)]

        return pairing_check(pairs)

    def verifiy_batch(self, proof: bytes, key: BatchBDFGVerifierKey) -> bool:

        transcript = self.new_transcript(proof)

        # read comitments
        commitments = [transcript.read_point() for _ in range(key.opening_size())]
        # get evaluation seed
        z = transcript.challenge()

        # read evaluations
        evals = [[transcript.read_scalar() for _ in range(opening.opening_size())] for opening in key.openings]

        # get combination base
        alpha = LinearCombination(transcript.challenge())

        # read first witness point
        W = transcript.read_point()

        # get linearisation point
        x = transcript.challenge()

        # read second witness point
        W_2 = transcript.read_point()

        # construct r_i(X) that is low degree equivalent of `f_i(X)`
        # where `f_i(set(z)) = r_i(set(z))`
        # then make a soft commitment to them `R_i = r_i(z) * G`

        r_i_x = [opening.low_degree_equivalent(evals, z) for opening, evals in zip(key.openings, evals)]
        # R_i = [self.G * r_x(z) for r_x in r_i_x]

        # pairing check:
        # `e(W', X) == e(z_x(z) * W' + L, G)``

        # first we reconstruct linearisation point `L`
        # ```
        # L = ∑ Z'_i(z) * (com(f_i(x)) - R_i)
        #   - com(h(x)) * Z(z)
        # ````

        # calculate `Z'_i(X) = Z _T(X) / Z_i(X)`
        z_i_x = [key.inverse_vanishing(i, z) for i in range(key.opening_size())]

        # `Z'_i(z) * (com(f_i(x)) - R_i)`
        linearisation_contibs = [(commitments[i] - self.G * r_i_x[i](x)) * z_i_x[i](x) for i in range(key.opening_size())]
        # combine linearisation
        L = alpha.combine_points(*linearisation_contibs)

        # `L = ∑ - com(h(x)) * Z(z)`
        z_x = key.vanishing(z)
        L = L - (W * z_x(x))

        # calcualate the escape from G2 term `z_x(z) * W'
        x_W_2 = W_2 * x

        # finally apply the pairing
        pairs = [(self.X_2, W_2), (self.n_G_2, x_W_2 + L)]

        return pairing_check(pairs)