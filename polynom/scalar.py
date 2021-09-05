from __future__ import annotations

one: Scalar = None
zero: Scalar = None
field_modulus = None
Scalar = None


def init_scalar_field(type):
    global Scalar
    Scalar = type
    global one
    one = Scalar(1)
    global zero
    zero = Scalar(0)
