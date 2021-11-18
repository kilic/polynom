from __future__ import annotations
from polynom.ecc import Scalar
from polynom.polynomial import Polynomial

# https://eprint.iacr.org/2020/081.pdf


def vanising_at(points: list[Scalar]):
    acc = Polynomial.one()
    for point in points:
        acc = acc * Polynomial.degree_one(point)
    return acc


class MultiBDFGCommon:

    def __init__(self, w: Scalar, shifts: list[int]):

        self.w = w
        self.shifts = shifts

    def opening_size(self) -> int:

        return len(self.shifts)

    def eval_points(self, z) -> list[Scalar]:

        return [z * (self.w**shift_val) for shift_val in self.shifts]

    def vanising(self, z) -> Polynomial:

        acc = Polynomial.one()

        for z in self.eval_points(z):
            acc = Polynomial.degree_one(z) * acc
        return acc


class BatchBDFGCommon:

    def eval_points(self, z) -> set:

        eval_points = set()
        for opening in self.openings:
            eval_points = eval_points | set(opening.eval_points(z))
        return eval_points

    def vanishing(self, z: Scalar):

        return vanising_at(self.eval_points(z))