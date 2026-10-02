"""The hand-worked answers. If these fail, trust nothing downstream."""
import numpy as np
import pytest

from limen import build_matrix, identifiability
from limen.toy import ALL, EXPECTED


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_expected_ranks(name):
    gene = ALL[name]()
    for read_length, want in EXPECTED[name].items():
        got = identifiability(gene, read_length).rank
        assert got == want, (
            f"{name} at read length {read_length}: rank {got}, expected {want}")


def test_the_central_claim():
    """Two distant splicing decisions are invisible to short reads.

    Four isoforms, two decisions about 2000 nt apart. Below the gap the matrix
    is rank deficient by exactly one: short reads measure two marginals and
    cannot recover the joint. Above the gap it is full rank."""
    gene = ALL["two_distant_skips"]()
    short = identifiability(gene, 150)
    long_ = identifiability(gene, 2300)

    assert short.n_iso == 4
    assert not short.identifiable
    assert short.deficiency == 1
    assert long_.identifiable
    assert long_.deficiency == 0


def test_distance_not_isoform_count_is_what_matters():
    """Same four isoforms, decisions 60 nt apart instead of 2000: resolved."""
    far = identifiability(ALL["two_distant_skips"](), 300)
    near = identifiability(ALL["two_adjacent_skips"](), 300)
    assert far.n_iso == near.n_iso == 4
    assert not far.identifiable
    assert near.identifiable


def test_different_populations_give_identical_data():
    """The rank deficiency is not an abstraction: two different mixtures
    produce identical observations below the gap, and differ above it."""
    gene = ALL["two_distant_skips"]()
    order = [t.id for t in gene.transcripts]

    a = np.zeros(4); a[order.index("ABMCD")] = 0.5; a[order.index("A_MD_")] = 0.5
    b = np.zeros(4); b[order.index("ABMD_")] = 0.5; b[order.index("A_MCD")] = 0.5

    M_short = build_matrix(gene, 150)
    M_long = build_matrix(gene, 2300)
    assert np.allclose(M_short @ a, M_short @ b)
    assert not np.allclose(M_long @ a, M_long @ b)


def test_full_length_read_is_always_identifiable():
    for name, build in ALL.items():
        assert identifiability(build(), 1_000_000).identifiable, name
