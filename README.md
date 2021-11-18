Experiments with polynomials and polynomial commitments

Includes various KZG based multi opening systems

* Original KZG
* Section 3 of [PLONK](https://eprint.iacr.org/2019/953.pdf)
* [Efficient polynomial commitment schemes for multiple points and polynomials](https://eprint.iacr.org/2020/081.pdf)

```python
def test_kzg_batch():

    n = 3
    KZG = kzg_setup(n)

    prover = KZG.prover_kzg()
    polys = [Polynomial.rand(1 << n) for _ in range(4)]
    proof = prover.create_proof_batch(polys)

    verifier = KZG.verifier_kzg()
    assert verifier.verify_batch(len(polys), proof)
```
