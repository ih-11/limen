"""The fast pattern finder must agree with brute force on every gene.

This is the test that matters. It caught three separate bugs during
development, each of which produced plausible but wrong numbers that the toy
genes alone did not expose."""
import random

import numpy as np
import pytest

from limen.identify import (Gene, Transcript, identifiability,
                            incidence_patterns, incidence_patterns_bruteforce)


def _random_gene(seed, n_exons=5, n_tx=4):
    rng = random.Random(seed)
    coords, pos = [], rng.randint(0, 50)
    for _ in range(n_exons):
        ln = rng.randint(20, 150)
        coords.append((pos, pos + ln))
        pos += ln + rng.randint(20, 400)
    txs, seen = [], set()
    for i in range(n_tx):
        keep = sorted(rng.sample(range(n_exons), rng.randint(1, n_exons)))
        ex = tuple(coords[k] for k in keep)
        if ex not in seen:
            seen.add(ex)
            txs.append(Transcript(f"t{i}", ex))
    return Gene(f"g{seed}", tuple(txs))


@pytest.mark.parametrize("seed", range(40))
def test_fast_matches_bruteforce_single_end(seed):
    gene = _random_gene(seed)
    for R in (20, 35, 60, 120, 300, 900):
        assert incidence_patterns(gene, R) == incidence_patterns_bruteforce(gene, R), \
            f"{gene.id} mismatch at read length {R}"


@pytest.mark.parametrize("seed", range(25))
def test_fast_matches_bruteforce_paired_end(seed):
    gene = _random_gene(1000 + seed)
    for R, F in ((20, 80), (35, 150), (50, 120), (30, 500)):
        assert incidence_patterns(gene, R, F) == \
               incidence_patterns_bruteforce(gene, R, F), \
            f"{gene.id} mismatch at read {R} fragment {F}"


@pytest.mark.parametrize("seed", range(20))
def test_rank_agrees_with_bruteforce(seed):
    gene = _random_gene(2000 + seed)
    for R in (30, 75, 200):
        fast = np.linalg.matrix_rank(
            np.array(sorted(incidence_patterns(gene, R)), dtype=float))
        slow = np.linalg.matrix_rank(
            np.array(sorted(incidence_patterns_bruteforce(gene, R)), dtype=float))
        assert int(fast) == int(slow)


def test_rank_never_exceeds_isoform_count():
    for seed in range(10):
        gene = _random_gene(3000 + seed)
        for R in (30, 100, 300):
            r = identifiability(gene, R)
            assert 0 <= r.rank <= r.n_iso


def test_single_isoform_is_trivially_identifiable():
    gene = Gene("g", (Transcript("only", ((0, 500),)),))
    assert identifiability(gene, 100).identifiable


def test_identical_structures_are_never_separable():
    """Two transcripts with the same exon chain can never be told apart."""
    ex = ((0, 100), (300, 400))
    gene = Gene("g", (Transcript("a", ex), Transcript("b", ex)))
    for R in (30, 100, 10_000):
        res = identifiability(gene, R)
        assert res.rank == 1 and res.deficiency == 1
