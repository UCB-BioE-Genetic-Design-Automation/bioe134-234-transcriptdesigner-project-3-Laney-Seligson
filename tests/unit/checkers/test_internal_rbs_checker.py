import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


@pytest.fixture
def checker():
    c = InternalRBSChecker()
    c.initiate()
    return c


def test_internal_rbs_detected(checker):
    seq = "GCTGCTAGGAGGAAAAATGGCTGCT"
    ok, bad_region = checker.run(seq)
    assert ok is False
    assert bad_region is not None


def test_no_internal_rbs(checker):
    seq = "GCTGCTGCTGCTGCTGCTGCTGCT"
    ok, bad_region = checker.run(seq)
    assert ok is True
    assert bad_region is None


def test_short_sequence_passes(checker):
    seq = "ATGGCT"
    ok, bad_region = checker.run(seq)
    assert ok is True
    assert bad_region is None