from polynom.ecc.bn254.scalar import Scalar as BN254_SCALAR
from polynom.ecc import init_scalar_field

init_scalar_field(BN254_SCALAR)

from polynom.ecc import init_ecc
from polynom.ecc.bn254.bn254 import BN254

init_ecc(BN254)

person = b"polynom test"
scalar_prefix = b"polynom test scalar"
challenge_prefix = b"polynom test challenge"
point_prefix = b"polynom test point"

from polynom.proof_system.transcript.transcript import TranscriptRead, TranscriptWrite
from polynom.proof_system.transcript.hasher import SHA256, Hasher


def hasher() -> Hasher:
    return SHA256(person, scalar_prefix, challenge_prefix, point_prefix)


def transcript_writer() -> TranscriptWrite:
    return TranscriptWrite(hasher())


def transcript_reader(message: bytes) -> TranscriptRead:
    return TranscriptRead(hasher(), message)
