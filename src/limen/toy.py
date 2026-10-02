"""
Toy genes with answers worked out by hand.

These are the project's ground truth.  If `identify.py` ever stops agreeing
with them, the bug is in the code and not in the biology, which is exactly
the check a simulation study needs before it is allowed to make any claim.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from .identify import Gene, Transcript


def _tx(name: str, parts: List[Tuple[str, int]], present: List[str]) -> Transcript:
    """Lay exons end to end in genomic space with 1 kb spacers between them."""
    exons, pos = [], 0
    for ex, ln in parts:
        if ex in present:
            exons.append((pos, pos + ln))
        pos += ln + 1000
    return Transcript(name, tuple(exons))


LAYOUT = [("A", 200), ("B", 100), ("M", 2000), ("C", 100), ("D", 200)]


def two_distant_skips() -> Gene:
    """
    A - [B] - M(2000) - [C] - D, with B and C independently skippable.

    Four isoforms; the two splicing decisions sit about 2000 nt apart in the
    mature message.  Any read shorter than that distance sees only the two
    marginals, what fraction contains B and what fraction contains C, and
    cannot recover which combinations actually existed.
    """
    combos = [("ABMCD", ["A", "B", "M", "C", "D"]),
              ("ABMD_", ["A", "B", "M", "D"]),
              ("A_MCD", ["A", "M", "C", "D"]),
              ("A_MD_", ["A", "M", "D"])]
    return Gene("two_distant_skips",
                tuple(_tx(n, LAYOUT, p) for n, p in combos))


def two_adjacent_skips() -> Gene:
    """Same four isoforms, but the skippable exons are 60 nt apart.

    The control for `two_distant_skips`: identical combinatorics, different
    geometry.  Short reads resolve this one, which shows the effect is about
    distance and not about isoform count."""
    layout = [("A", 200), ("B", 100), ("M", 60), ("C", 100), ("D", 200)]
    combos = [("ABMCD", ["A", "B", "M", "C", "D"]),
              ("ABMD_", ["A", "B", "M", "D"]),
              ("A_MCD", ["A", "M", "C", "D"]),
              ("A_MD_", ["A", "M", "D"])]
    return Gene("two_adjacent_skips",
                tuple(_tx(n, layout, p) for n, p in combos))


def single_cassette() -> Gene:
    """One skipped exon, two isoforms. Identifiable at any usable read length."""
    layout = [("A", 200), ("B", 100), ("C", 200)]
    return Gene("single_cassette",
                (_tx("inc", layout, ["A", "B", "C"]),
                 _tx("skp", layout, ["A", "C"])))


def nested_ends() -> Gene:
    """
    Two isoforms sharing every junction, differing only by a longer 3' exon.

    Structurally identifiable at any read length, because the extension is
    unique to the long isoform.  Included as a contrast: identifiable is not
    the same as well determined.  At finite depth very few reads land in the
    distinguishing region, which is the regime where quantifiers start
    inventing answers.
    """
    return Gene("nested_ends",
                (Transcript("long",  ((0, 200), (1200, 1900))),
                 Transcript("short", ((0, 200), (1200, 1400)))))


#: gene name -> {read length: expected rank}
EXPECTED: Dict[str, Dict[int, int]] = {
    # B and C sit 2000 nt apart; a read must run from inside B to inside C,
    # so it needs 1 + 2000 + 1 = 2002 nt before the fourth dimension appears.
    "two_distant_skips":  {60: 3, 150: 3, 500: 3, 2000: 3, 2100: 4, 5000: 4},
    # identical combinatorics, 60 nt apart: 62 nt suffices.
    "two_adjacent_skips": {60: 3, 100: 4, 300: 4, 5000: 4},
    "single_cassette":    {60: 2, 100: 2, 300: 2},
    "nested_ends":        {60: 2, 100: 2, 300: 2},
}

ALL = {
    "two_distant_skips": two_distant_skips,
    "two_adjacent_skips": two_adjacent_skips,
    "single_cassette": single_cassette,
    "nested_ends": nested_ends,
}
