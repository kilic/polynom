from __future__ import annotations, barry_as_FLUFL
from polynom.commitment.kzg_base import KZGProverBase
from polynom.commitment.bdfg.common import BatchBDFGCommon, MultiBDFGCommon, vanising_at
from polynom.ecc import Scalar, zero
from polynom.polynomial import Polynomial, lagrange_interpolation
from polynom.domain import Domain
from polynom.lc import LinearCombination


class MultiBDFGProverKey(MultiBDFGCommon):

    def __init__(self, domain: Domain, poly: Polynomial, shifts: list[int]):

        super().__init__(domain.w(), shifts)
        self.domain = domain
        self.poly = poly

    def evaluate(self, z) -> list[Scalar]:

        return [self.poly(e) for e in self.eval_points(z)]

    def low_degree_equivalent(self, z: Scalar) -> Polynomial:

        eval_points = self.eval_points(z)
        evals = self.evaluate(z)
        points = [(x, y) for x, y in zip(eval_points, evals)]
        return lagrange_interpolation(points)

    def linearision_polynomial(self, z: Scalar) -> Polynomial:

        z_x = self.vanising(z)
        z_eval = z_x(z)
        r_x = self.low_degree_equivalent(z)
        r_eval = r_x(z)
        h_x = self.quotient_polynomial(z)
        l_x = self.poly - r_eval - h_x * z_eval
        # sanity check
        assert l_x(z) == zero
        return l_x

    def quotient_polynomial(self, z: Scalar) -> Polynomial:

        z_x = self.vanising(z)
        r_x = self.low_degree_equivalent(z)
        q_x = self.domain.div(self.poly - r_x, z_x)
        # sanity check
        assert q_x.degree() == self.poly.degree() - z_x.degree()
        return q_x

    def linearized_quotient_polynomial(self, z) -> Polynomial:

        l_x = self.linearision_polynomial(z)
        u_x = self.domain.div(l_x, Polynomial.degree_one(z))
        # sanity check
        assert u_x.degree() + 1 == l_x.degree()
        return u_x


class BatchBDFGProverKey(BatchBDFGCommon):

    def __init__(self, domain: Domain, openings: list[MultiBDFGProverKey]):
        self.openings = openings
        self.domain = domain

    def multi_open_size(self) -> int:
        return len(self.openings)

    def polynomials(self) -> list[Polynomial]:

        return [opening.poly for opening in self.openings]

    def evaluate(self, z: Scalar) -> list[list[Scalar]]:

        return [opening.evaluate(z) for opening in self.openings]

    def inverse_vanishing(self, multi_open_index: int, z: Scalar) -> Polynomial:

        acc = Polynomial.one()
        # TODO: double check the order
        inv_eval_points = self.eval_points(z) ^ set(self.openings[multi_open_index].eval_points(z))
        return vanising_at(inv_eval_points)

    # # TODO: combine all here
    # def combine_linearisation_contribs(self, multi_open_index: int, z: Scalar) -> Polynomial:
    #     # [((self.openings[i].poly - self.openings[i].low_degree_equivalent(z)) * self.inverse_vanishing(i, z)) for i in range(len(self.openings))]

    def linearisation_contrib(self, multi_open_index: int, z: Scalar) -> Polynomial:

        z_i_not_x = self.inverse_vanishing(multi_open_index, z)
        f_x = self.openings[multi_open_index].poly
        r_x = self.openings[multi_open_index].low_degree_equivalent(z)
        return (f_x - r_x(z)) * z_i_not_x(z)

    def quotient_polynomial(self, alpha: LinearCombination, z: Scalar):

        quotient_contribs = [opening.quotient_polynomial(z) for opening in self.openings]
        return alpha.combine_poly(*quotient_contribs)

    def linearized_quotient_polynomial(self, alpha: LinearCombination, z: Scalar):

        linearisation_contribs = [self.linearisation_contrib(i, z) for i in range(len(self.openings))]
        q_x = self.quotient_polynomial(alpha, z)
        z_x = self.vanishing(z)
        z_eval = z_x(z)
        l_x = alpha.combine_poly(*linearisation_contribs) - q_x * z_eval
        # sanity check
        assert l_x(z) == zero

        return self.domain.div(l_x, Polynomial.degree_one(z))


class BDFGProver(KZGProverBase):

    def new_multi_key(self, poly: Polynomial, shifts: list[int]) -> MultiBDFGProverKey:
        return MultiBDFGProverKey(self.domain, poly, shifts)

    def new_batch_key(self, openings: list[MultiBDFGProverKey]) -> BatchBDFGProverKey:
        return BatchBDFGProverKey(self.domain, openings)

    def create_proof_single(self, multi_point: MultiBDFGProverKey) -> bytes:

        transcript = self.new_transcript()

        # commit to the polynomial f(X) and write it to the trasctipy
        transcript.write_point(self.commit(multi_point.poly))

        # get evaluation point
        z = transcript.challenge()

        # evaluate polynomial at shifted values of z
        # [f(z_0), f(z_1), ...]
        # then write them to the transcript
        [transcript.write_scalar(eval) for eval in multi_point.evaluate(z)]

        # calculate first quotient h(X)
        # commit to the quotient W = com(h(X))
        # write it to the transcript
        transcript.write_point(self.commit(multi_point.quotient_polynomial(z)))
        # same goes for second quotient h'(Xf)
        transcript.write_point(self.commit(multi_point.linearized_quotient_polynomial(z)))
        return transcript.get_message()

    def create_proof_batch(self, batch: BatchBDFGProverKey) -> bytes:

        transcript = self.new_transcript()

        # commit to polynomials f_i(X) and write commitments to the transcript
        [transcript.write_point(c) for c in self.c(*batch.polynomials())]
        # get evaluation point
        z = transcript.challenge()

        # evaluate polynomials in multiple points
        # [f_0(u_0), f_0(u_1), ...]
        # [f_1(v_0), f_1(v_1), ...]
        evals = batch.evaluate(z)
        # write evaluations to the transcript
        [[transcript.write_scalar(eval) for eval in evals_i] for evals_i in evals]

        # get combination base
        alpha = LinearCombination(transcript.challenge())

        # calculate first quotient h(X)
        h_x = batch.quotient_polynomial(alpha, z)
        # commit to the first quotient and write it to the transcript
        transcript.write_point(self.commit(h_x))

        # calculate second quotient
        h2_x = batch.linearized_quotient_polynomial(alpha, z)
        # commit to the scond quotient and write it to the transcript
        transcript.write_point(self.commit(h2_x))

        # return proof
        return transcript.get_message()