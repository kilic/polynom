from polynom.commitment import gwc
from polynom.commitment.bdfg.prover import BatchBDFGProverKey, MultiBDFGProverKey
from polynom.commitment.gwc import GWCKey
from polynom.commitment.kzg_setup import KZGSetup
from polynom.polynomial import Polynomial
from polynom.ecc import Scalar
from polynom.domain import Domain
from polynom.ecc.bn254.domain import domain_config
from . import hasher


def kzg_setup(n: int) -> KZGSetup:
    domain = Domain(domain_config(n))
    return KZGSetup.new(domain, hasher())


def test_kzg_single():

    n = 3
    KZG = kzg_setup(n)

    prover = KZG.prover_kzg()
    p_x = Polynomial.rand(1 << n)
    proof = prover.create_proof(p_x)

    verifier = KZG.verifier_kzg()
    assert verifier.verify(proof)


def test_kzg_batch():

    n = 3
    KZG = kzg_setup(n)

    prover = KZG.prover_kzg()
    polys = [Polynomial.rand(1 << n) for _ in range(4)]
    proof = prover.create_proof_batch(polys)

    verifier = KZG.verifier_kzg()
    assert verifier.verify_batch(len(polys), proof)


def test_kzg_multi_gwc():

    n = 3
    KZG = kzg_setup(n)

    polys = [Polynomial.rand(1 << n) for _ in range(4)]

    shift_map = {0: [1, 0], 1: [0, 1, 3, 2], 2: [3], 19: [3, 2, 1, 0]}
    prover = KZG.prover_gwc()
    context = GWCKey(shift_map)
    proof = prover.create_proof(polys, context)

    verifier = KZG.verifier_gwc()
    assert verifier.verify(context, proof)


def test_bdfg_precomputation():

    n = 6
    KZG = kzg_setup(n)
    domain = KZG.domain
    zero = Scalar(0)

    p_0_x = Polynomial.rand(1 << n)
    shifts_0 = [0, 1, 2, 3, 10, 11]
    multi_0 = MultiBDFGProverKey(domain, p_0_x, shifts_0)
    z = Scalar.rand()

    r_x = multi_0.low_degree_equivalent(z)
    for z in multi_0.eval_points(z):
        assert p_0_x(z) == r_x(z)
    assert r_x.degree() == len(shifts_0) - 1

    z_x = multi_0.vanising(z)
    assert z_x.degree() == len(shifts_0)
    for u in multi_0.eval_points(z):
        assert z_x(u) == zero

    q_x = multi_0.quotient_polynomial(z)
    assert q_x.degree() == p_0_x.degree() - len(shifts_0)

    l_x = multi_0.linearized_quotient_polynomial(z)
    assert l_x.degree() == p_0_x.degree() - 1

    p_1_x = Polynomial.rand(1 << n)
    shifts_1 = [1, 2, 4, 5]
    multi_1 = MultiBDFGProverKey(domain, p_1_x, shifts_1)

    batch = BatchBDFGProverKey(domain, [multi_0, multi_1])

    eval_points = batch.eval_points(z)
    assert len(eval_points) == len(set(shifts_0) | set(shifts_1))

    eval_points = set(multi_1.eval_points(z)) - set(multi_0.eval_points(z))
    z_inv_0_x = batch.inverse_vanishing(0, z)
    for e in eval_points:
        z_inv_0_x(e) == zero

    eval_points = set(multi_0.eval_points(z)) - set(multi_1.eval_points(z))
    z_inv_1_x = batch.inverse_vanishing(1, z)
    for e in eval_points:
        z_inv_1_x(e) == zero

    eval_points = batch.eval_points(z)
    len(eval_points) == len(set(multi_0.eval_points(z)) | set(multi_1.eval_points(z)))
    z_x = batch.vanishing(z)
    assert z_x.degree() == len(eval_points)
    for u in eval_points:
        assert z_x(u) == zero


def test_bdfg_single():

    n = 6
    KZG = kzg_setup(n)
    prover, verifier = KZG.prover_bdfg(), KZG.verifier_bdfg()

    p_0_x = Polynomial.rand(1 << n)
    shifts_0 = [1, 2, 3, 10, 11]
    prover_key = prover.new_multi_key(p_0_x, shifts_0)
    proof = prover.create_proof_single(prover_key)

    verifier_key = verifier.new_multi_key(shifts_0)
    assert verifier.verify_single(proof, verifier_key)


def test_bdfg_batch():

    n = 6
    KZG = kzg_setup(n)
    prover, verifier = KZG.prover_bdfg(), KZG.verifier_bdfg()

    p_0_x = Polynomial.rand(1 << n)
    shifts_0 = [1, 2, 3]
    key_multi_0 = prover.new_multi_key(p_0_x, shifts_0)

    p_1_x = Polynomial.rand(1 << n)
    shifts_1 = [1, 2]
    key_multi_1 = prover.new_multi_key(p_1_x, shifts_1)
    batch = prover.new_batch_key([key_multi_0, key_multi_1])
    proof = prover.create_proof_batch(batch)

    key_multi_0, key_multi_1 = verifier.new_multi_key(shifts_0), verifier.new_multi_key(shifts_1)
    key = verifier.new_batch_key([key_multi_0, key_multi_1])
    assert verifier.verifiy_batch(proof, key)


def test_bdfg_batch_mal():

    n = 6
    KZG = kzg_setup(n)
    prover, verifier = KZG.prover_bdfg(), KZG.verifier_bdfg()

    p_0_x = Polynomial.rand(1 << n)
    shifts_0 = [1, 2]
    key_multi_0 = prover.new_multi_key(p_0_x, shifts_0)

    p_1_x = Polynomial.rand(1 << n)

    shifts_1 = [0, 3]  # this one makes malicious case work!

    key_multi_1 = prover.new_multi_key(p_1_x, shifts_1)
    batch = prover.new_batch_key([key_multi_0, key_multi_1])
    proof = prover.create_proof_batch_mal(batch)

    key_multi_0, key_multi_1 = verifier.new_multi_key(shifts_0), verifier.new_multi_key(shifts_1)
    key = verifier.new_batch_key([key_multi_0, key_multi_1])
    assert verifier.verifiy_batch(proof, key)
