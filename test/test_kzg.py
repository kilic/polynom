from polynom.commitment.kzg import BatchGWC, KZGSetup, KZGProver, KZGVerifier
from polynom.polynomial import Polynomial
from polynom.ecc import Scalar
from polynom.lc import LinearCombination
from polynom.domain import Domain
from polynom.ecc.bn254.domain import domain_config
from . import hasher


def kzg_setup(n: int) -> KZGSetup:
    domain = Domain(domain_config(n))
    return KZGSetup.new(domain, hasher())


def test_kzg_single_poly():

    n = 3
    KZG = kzg_setup(n)

    prover = KZG.prover()
    p_x = Polynomial.rand(1 << n)
    proof = prover.create_proof(p_x)

    verifier = KZG.verifier()
    assert verifier.verify(proof)


def test_kzg_batch_single_open():

    n = 3
    KZG = kzg_setup(n)

    N = 4
    prover = KZG.prover()
    polys = [Polynomial.rand(4) for _ in range(N)]
    proof = prover.create_proof_single_point_open_batch(polys)

    verifier = KZG.verifier()
    assert verifier.verify_single_point_open_batch(len(polys), proof)


def test_kzg_batch_multi_open():

    n = 3
    KZG = kzg_setup(n)

    N = 4
    polys = [Polynomial.rand(4) for _ in range(N)]

    shift_map = {0: [1, 0], 1: [0, 1, 3, 2], 2: [3], 19: [3, 2, 1, 0]}
    prover = KZG.prover()
    context = BatchGWC(shift_map)
    proof = prover.create_proof_multi_point_open_batch_gwc(polys, context)

    verifier = KZG.verifier()
    assert verifier.verify_multi_point_open_batch_gwc(context, proof)
