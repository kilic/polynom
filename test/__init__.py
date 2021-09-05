from polynom.ecc.bn254.scalar import Scalar as BN254_SCALAR
from polynom.scalar import init_scalar_field

init_scalar_field(BN254_SCALAR)

from polynom.ecc.curve import init_ecc
from polynom.ecc.bn254.bn254 import BN254

init_ecc(BN254)
