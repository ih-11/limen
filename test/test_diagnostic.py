"""The fast conditioning path must equal brute force on every gene."""
import random

import pytest

from limen.identify import (Gene, Transcript, diagnostic_fraction,
                            signature_owners, signature_owners_bruteforce)


def _diagnostic_bruteforce(gene, R, F=None):
    K = gene.n_iso
    priv, tot = [0] * K, [0] * K
    for _s, who in signature_owners_bruteforce(gene, R, F).items():
        solo = len(who) == 1
        for j, n in who.items():
            tot[j] += n
            if solo:
                priv[j] += n
    return list(zip(priv, tot))


def _gene(seed, n_ex=6, n_tx=5):
    rng = random.Random(seed)
    coords, pos = [], rng.randint(0, 50)
    for _ in range(n_ex):
        ln = rng.randint(15, 400)
        coords.append((pos, pos + ln))
        pos += ln + rng.randint(1, 600)
    txs, seen = [], set()
    for i in range(n_tx):
        keep = sorted(rng.sample(range(n_ex), rng.randint(1, n_ex)))
        ex = tuple(coords[k] for k in keep)
        if ex not in seen:
            seen.add(ex)
            txs.append(Transcript(f"t{i}", ex))
    return Gene(f"g{seed}", tuple(txs))


@pytest.mark.parametrize("seed", range(30))
def test_single_end(seed):
    g = _gene(seed)
    for R in (10, 25, 60, 150, 400, 1500, 9999):
        assert diagnostic_fraction(g, R) == _diagnostic_bruteforce(g, R), f"R={R}"


@pytest.mark.parametrize("seed", range(20))
def test_paired_end(seed):
    g = _gene(1000 + seed)
    for R, F in ((25, 100), (50, 200), (75, 400), (150, 150), (300, 1000)):
        assert diagnostic_fraction(g, R, F) == _diagnostic_bruteforce(g, R, F)


@pytest.mark.parametrize("shape", [
    (((0, 5), (10, 15)), ((0, 5), (20, 25))),
    (((0, 100),), ((0, 100),)),
    (((0, 50), (50, 100)), ((0, 100),)),
    (((0, 1), (5, 6), (10, 11)), ((0, 1), (10, 11))),
])
def test_edge_cases(shape):
    g = Gene("e", tuple(Transcript(f"t{i}", ex) for i, ex in enumerate(shape)))
    for R in (1, 2, 3, 7, 50, 500):
        assert diagnostic_fraction(g, R) == _diagnostic_bruteforce(g, R)


def test_identical_transcripts_have_no_diagnostic_evidence():
    ex = ((0, 100), (300, 400))
    g = Gene("g", (Transcript("a", ex), Transcript("b", ex)))
    for R in (20, 150, 5000):
        assert all(p == 0 for p, _ in diagnostic_fraction(g, R))
