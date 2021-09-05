from polynom.commitment.kzg import KZGProof, kzg_setup
from polynom.polynomial import Polynomial
from polynom.scalar import Scalar
from polynom.lc import LinearCombination
from polynom.domain import Domain
from polynom.ecc.bn254.domain import domain_config


def test_commitment():
    n = 3
    domain = Domain(domain_config(n))
    prover, _ = kzg_setup(1 << n, domain)
    v = [Polynomial.rand(1 << n) for _ in range(10)]
    vi0 = prover.c(*v)
    vi1 = [prover.commit(u) for u in v]
    for v0, v1 in zip(vi0, vi1):
        assert v0 == v1


def test_kzg():

    n = 3
    domain = Domain(domain_config(n))
    prover, verifier = kzg_setup(1 << n, domain)

    p_x = Polynomial.rand(1 << n)

    z = Scalar.rand()
    C = prover.commit(p_x)
    e = p_x(z)
    W = prover.witness(z, p_x)
    proof = KZGProof(z, W, C, e)
    assert verifier.verify(proof)


def test_kzg_batch():

    n = 3
    domain = Domain(domain_config(n))
    prover, verifier = kzg_setup(1 << n, domain)

    p0_x_z0, p1_x_z0, p2_x_z0 = Polynomial.rand(4), Polynomial.rand(4), Polynomial.rand(4)
    p0_x_z1, p1_x_z1, p2_x_z1 = Polynomial.rand(4), Polynomial.rand(4), Polynomial.rand(4)
    p0_x_z2, p1_x_z2, p2_x_z2 = Polynomial.rand(4), Polynomial.rand(4), Polynomial.rand(4)
    z0, z1, z2 = Scalar.rand(), Scalar.rand(), Scalar.rand()

    ch0 = LinearCombination(Scalar.rand())
    ch1 = LinearCombination(Scalar.rand())
    ch2 = LinearCombination(Scalar.rand())

    input = [
        (ch0, z0, [p0_x_z0, p1_x_z0, p2_x_z0]),
        (ch1, z1, [p0_x_z1, p1_x_z1, p2_x_z1]),
        (ch2, z2, [p0_x_z2, p1_x_z2, p2_x_z2]),
    ]

    P0 = prover.c(p0_x_z0, p1_x_z0, p2_x_z0)
    P1 = prover.c(p0_x_z1, p1_x_z1, p2_x_z1)
    P2 = prover.c(p0_x_z2, p1_x_z2, p2_x_z2)

    e0 = [p0_x_z0(z0), p1_x_z0(z0), p2_x_z0(z0)]
    e1 = [p0_x_z1(z1), p1_x_z1(z1), p2_x_z1(z1)]
    e2 = [p0_x_z2(z2), p1_x_z2(z2), p2_x_z2(z2)]

    W = prover.witness_batch(input)

    proof = [
        KZGProof(z0, W[0], ch0.combine_points(*P0), ch0.combine_fr(*e0)),
        KZGProof(z1, W[1], ch1.combine_points(*P1), ch1.combine_fr(*e1)),
        KZGProof(z2, W[2], ch2.combine_points(*P2), ch2.combine_fr(*e2)),
    ]

    assert verifier.verify_batch(LinearCombination(Scalar.rand()), proof)
