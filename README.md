Experiments with polynomials and polynomial commitments.

```python
def test():

    n = 3
    domain = Domain(domain_config(n))
    prover, verifier = kzg_setup(1 << n, domain)

    p_x = Polynomial.rand(1 << n)

    assert domain.evaluate(domain.interpolate(p_x)) == p_x
    assert domain.interpolate(domain.evaluate(p_x)) == p_x

    z = Scalar.rand()
    C = prover.commit(p_x)
    e = p_x(z)
    W = prover.witness(z, p_x)
    proof = KZGProof(z, W, C, e)
    assert verifier.verify(proof)
```
