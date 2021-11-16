from polynom.ecc import Point, Scalar, CURVE
from polynom.proof_system.transcript.transcript import Transcript, TranscriptRead, TranscriptWrite
from polynom.proof_system.transcript.hasher import SHA256
from . import hasher


def test_conversion():
    p0 = Point.rand()
    uncompressed = p0.to_uncompressed()
    p1 = Point.from_uncompressed(uncompressed)
    assert p0 == p1


def test_transcript():

    p0 = Point.rand()
    p1 = Point.rand()
    p2 = Point.rand()
    p_common = Point.rand()
    s0 = Scalar.rand()
    s1 = Scalar.rand()
    s2 = Scalar.rand()
    s_common = Scalar.rand()

    t = TranscriptWrite(hasher())
    t.write_point(p0)
    ch_0 = t.challenge()
    t.write_scalar(s0)
    ch_1 = t.challenge()
    t.write_point(p1)
    ch_2 = t.challenge()
    t.write_point(p2)
    ch_3 = t.challenge()
    t.write_scalar(s1)
    ch_4 = t.challenge()
    t.write_scalar(s2)
    ch_5 = t.challenge()
    t.write_point_to_state(p_common)
    ch_6 = t.challenge()
    t.write_scalar_to_state(s_common)
    ch_7 = t.challenge()

    message = t.get_message()

    t = TranscriptRead(hasher(), message)
    p0_read = t.read_point()
    ch_0_read = t.challenge()
    s0_read = t.read_scalar()
    ch_1_read = t.challenge()
    p1_read = t.read_point()
    ch_2_read = t.challenge()
    p2_read = t.read_point()
    ch_3_read = t.challenge()
    s1_read = t.read_scalar()
    ch_4_read = t.challenge()
    s2_read = t.read_scalar()
    ch_5_read = t.challenge()
    t.write_point_to_state(p_common)
    ch_6_read = t.challenge()
    t.write_scalar_to_state(s_common)
    ch_7_read = t.challenge()

    assert p0_read == p0
    assert p1_read == p1
    assert p2_read == p2
    assert s0_read == s0
    assert s1_read == s1
    assert s2_read == s2

    assert ch_0 == ch_0_read
    assert ch_1 == ch_1_read
    assert ch_2 == ch_2_read
    assert ch_3 == ch_3_read
    assert ch_4 == ch_4_read
    assert ch_5 == ch_5_read
    assert ch_6 == ch_6_read
    assert ch_7 == ch_7_read
