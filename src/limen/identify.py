"""
Level 0: structural identifiability of transcript isoforms.

A sequencing read observes a contiguous stretch of one mature transcript.
Once aligned, all that remains visible is the set of genomic blocks it
covers.  Two reads from different transcripts are indistinguishable exactly
when those block sets coincide.

For a gene with K transcripts, every possible read gives one row of a 0/1
incidence matrix M: a 1 for each transcript that could have produced that
read.  Entries are 0/1 because a block set fixes the genomic start, which
fixes the position within any transcript able to produce it.  Then

    rank(M) == K   ->  abundances are uniquely determined: IDENTIFIABLE
    rank(M) <  K   ->  genuinely different transcript mixtures produce
                       byte-identical data, and no algorithm, no depth and
                       no replication can separate them.

K - rank(M) is the number of invisible dimensions.

This is the infinite-depth, error-free ceiling.  Real pipelines land below
it; that gap is what Level 1 measures.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

Block = Tuple[int, int]            # genomic interval, half-open [start, end)
Signature = Tuple                  # blocks (single-end) or a pair of those
Pattern = Tuple[int, ...]          # one row of M

__all__ = ["Transcript", "Gene", "GeneResult", "read_signature",
           "incidence_patterns", "incidence_patterns_bruteforce",
           "build_matrix", "identifiability", "sweep"]


# --------------------------------------------------------------------------
# data model
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Transcript:
    """Exons are genomic, half-open, ascending and non-overlapping."""
    id: str
    exons: Tuple[Block, ...]

    @property
    def length(self) -> int:
        return sum(e - s for s, e in self.exons)

    @property
    def boundaries(self) -> Tuple[int, ...]:
        """Internal exon ends, in TRANSCRIPT coordinates."""
        out, off = [], 0
        for s, e in self.exons[:-1]:
            off += e - s
            out.append(off)
        return tuple(out)

    def to_genomic(self, a: int, b: int) -> Tuple[Block, ...]:
        """Map transcript interval [a, b) onto genomic blocks."""
        blocks: List[Block] = []
        off = 0
        for s, e in self.exons:
            ln = e - s
            lo, hi = max(a, off), min(b, off + ln)
            if hi > lo:
                blocks.append((s + (lo - off), s + (hi - off)))
            off += ln
            if off >= b:
                break
        return tuple(blocks)

    def to_transcript(self, g: int) -> Optional[int]:
        """Transcript coordinate of genomic position g, or None if not exonic."""
        off = 0
        for s, e in self.exons:
            if s <= g < e:
                return off + (g - s)
            off += e - s
        return None


@dataclass(frozen=True)
class Gene:
    id: str
    transcripts: Tuple[Transcript, ...]
    seqid: str = ""
    strand: str = "."

    @property
    def n_iso(self) -> int:
        return len(self.transcripts)


@dataclass(frozen=True)
class GeneResult:
    gene_id: str
    n_iso: int
    rank: int
    read_length: int
    fragment_length: Optional[int]
    n_patterns: int

    @property
    def identifiable(self) -> bool:
        return self.rank == self.n_iso

    @property
    def deficiency(self) -> int:
        return self.n_iso - self.rank


# --------------------------------------------------------------------------
# what a single read looks like
# --------------------------------------------------------------------------

def _windows(x: int, R: int, F: Optional[int]) -> Tuple[Tuple[int, int], ...]:
    """Transcript-coordinate windows actually sequenced, for fragment start x."""
    if F is None:
        return ((x, x + R),)
    if F <= R:                                   # mates cover the whole fragment
        return ((x, x + F),)
    return ((x, x + R), (x + F - R, x + F))


def read_signature(tx: Transcript, x: int, R: int,
                   F: Optional[int] = None) -> Signature:
    """Genomic footprint of a read whose fragment starts at transcript coord x."""
    parts = tuple(tx.to_genomic(a, b) for a, b in _windows(x, R, F))
    return parts[0] if len(parts) == 1 else parts


def _whole(tx: Transcript, F: Optional[int]) -> Signature:
    w = tx.to_genomic(0, tx.length)
    return w if F is None else (w, w)


# --------------------------------------------------------------------------
# the incidence patterns
# --------------------------------------------------------------------------

def _offsets(R: int, F: Optional[int]) -> Tuple[int, ...]:
    return (0, R) if F is None else (0, R, F - R, F)


def incidence_patterns(gene: Gene, R: int,
                       F: Optional[int] = None) -> Set[Pattern]:
    """
    Distinct rows of M, found without enumerating every read.

    Which transcripts share a read's footprint can only change where some
    transcript's read window crosses an exon boundary or runs off its end.
    Between those genomic breakpoints the grouping is constant, so sampling
    one position per interval is exact.  Cost is O(exons x transcripts)
    rather than O(transcript length).
    """
    txs = gene.transcripts
    span = R if F is None else F

    cuts: Set[int] = set()
    for t in txs:
        if not t.exons:
            continue
        L = t.length
        if span >= L:
            cuts.add(t.exons[0][0])
            continue
        max_x = L - span
        xs = {0, max_x}
        for b in t.boundaries:
            for o in _offsets(R, F):
                x = b - o
                if 0 <= x <= max_x:
                    xs.add(x)
                if 0 <= x - 1 <= max_x:
                    xs.add(x - 1)
        for x in xs:
            blk = t.to_genomic(x, x + 1)
            if blk:
                cuts.add(blk[0][0])

    patterns: Set[Pattern] = set()
    for g in sorted(cuts):
        groups: Dict[Signature, List[int]] = {}
        for j, t in enumerate(txs):
            if not t.exons:
                continue
            L = t.length
            span_t = R if F is None else F
            if span_t >= L:
                if g == t.exons[0][0]:
                    groups.setdefault(_whole(t, F), []).append(j)
                continue
            x = t.to_transcript(g)
            if x is None or x > L - span_t:
                continue
            groups.setdefault(read_signature(t, x, R, F), []).append(j)
        for js in groups.values():
            row = [0] * len(txs)
            for j in js:
                row[j] = 1
            patterns.add(tuple(row))
    return patterns


def incidence_patterns_bruteforce(gene: Gene, R: int,
                                  F: Optional[int] = None) -> Set[Pattern]:
    """Reference implementation: every signature of every transcript. Slow."""
    txs = gene.transcripts
    owners: Dict[Signature, List[int]] = {}
    for j, t in enumerate(txs):
        if not t.exons:
            continue
        L = t.length
        span = R if F is None else F
        if span >= L:
            owners.setdefault(_whole(t, F), []).append(j)
            continue
        for x in range(0, L - span + 1):
            owners.setdefault(read_signature(t, x, R, F), []).append(j)
    patterns: Set[Pattern] = set()
    for js in owners.values():
        row = [0] * len(txs)
        for j in js:
            row[j] = 1
        patterns.add(tuple(row))
    return patterns


# --------------------------------------------------------------------------
# rank
# --------------------------------------------------------------------------

def build_matrix(gene: Gene, R: int, F: Optional[int] = None) -> np.ndarray:
    """The 0/1 incidence matrix M, duplicate rows already collapsed."""
    pats = incidence_patterns(gene, R, F)
    if not pats:
        return np.zeros((0, gene.n_iso), dtype=float)
    return np.array(sorted(pats), dtype=float)


def identifiability(gene: Gene, R: int,
                    F: Optional[int] = None) -> GeneResult:
    """Is this gene structurally identifiable at this read length?"""
    if gene.n_iso == 0:
        return GeneResult(gene.id, 0, 0, R, F, 0)
    M = build_matrix(gene, R, F)
    rank = int(np.linalg.matrix_rank(M)) if M.size else 0
    return GeneResult(gene.id, gene.n_iso, rank, R, F, M.shape[0])


def sweep(genes: Iterable[Gene], read_lengths: Sequence[int],
          fragment_lengths: Optional[Sequence[Optional[int]]] = None
          ) -> List[GeneResult]:
    """identifiability over a grid of genes x read lengths x fragment lengths."""
    frags: Sequence[Optional[int]] = (list(fragment_lengths)
                                      if fragment_lengths is not None else [None])
    out: List[GeneResult] = []
    for g in genes:
        for R in read_lengths:
            for F in frags:
                if F is not None and F < R:
                    continue
                out.append(identifiability(g, R, F))
    return out


# --------------------------------------------------------------------------
# conditioning: not "is it identifiable" but "how much evidence is there"
# --------------------------------------------------------------------------

def signature_owners(gene: Gene, R: int,
                     F: Optional[int] = None) -> Dict[Signature, Dict[int, int]]:
    """
    {signature: {transcript index: number of start positions producing it}}.

    Rank asks whether a distinguishing read exists.  This asks how many there
    are, which is the question that survives contact with finite depth.
    """
    txs = gene.transcripts
    span = R if F is None else F

    cuts: Set[int] = set()
    for t in txs:
        if not t.exons:
            continue
        for s, e in t.exons:                   # every exon edge bounds a segment
            cuts.add(s)
            cuts.add(e)
        L = t.length
        if span >= L:
            continue
        max_x = L - span
        for b in t.boundaries:                 # and every window-edge crossing
            for o in _offsets(R, F):
                for x in (b - o - 1, b - o, b - o + 1):
                    if 0 <= x <= max_x:
                        blk = t.to_genomic(x, x + 1)
                        if blk:
                            cuts.add(blk[0][0])
        for x in (0, max_x, max_x + 1):        # including the first invalid start
            blk = t.to_genomic(x, x + 1)
            if blk:
                cuts.add(blk[0][0])
    ordered = sorted(cuts)

    owners: Dict[Signature, Dict[int, int]] = {}
    for i, g in enumerate(ordered):
        stop = ordered[i + 1] if i + 1 < len(ordered) else g + 1
        w = stop - g
        if w <= 0:
            continue
        for j, t in enumerate(txs):
            if not t.exons:
                continue
            L = t.length
            if span >= L:
                continue
            x = t.to_transcript(g)
            if x is None or x > L - span:
                continue
            sig = read_signature(t, x, R, F)
            owners.setdefault(sig, {})
            owners[sig][j] = owners[sig].get(j, 0) + w

    for j, t in enumerate(txs):                # reads longer than the molecule
        if t.exons and span >= t.length:
            sig = _whole(t, F)
            owners.setdefault(sig, {})
            owners[sig][j] = owners[sig].get(j, 0) + 1
    return owners


def signature_owners_bruteforce(gene: Gene, R: int,
                                F: Optional[int] = None) -> Dict[Signature, Dict[int, int]]:
    """Reference implementation.  Slow; used only to validate the fast path."""
    txs = gene.transcripts
    span = R if F is None else F
    owners: Dict[Signature, Dict[int, int]] = {}
    for j, t in enumerate(txs):
        if not t.exons:
            continue
        if span >= t.length:
            sig = _whole(t, F)
            owners.setdefault(sig, {})
            owners[sig][j] = owners[sig].get(j, 0) + 1
            continue
        for x in range(0, t.length - span + 1):
            sig = read_signature(t, x, R, F)
            owners.setdefault(sig, {})
            owners[sig][j] = owners[sig].get(j, 0) + 1
    return owners


def private_fraction(gene: Gene, R: int,
                     F: Optional[int] = None) -> List[Tuple[int, int]]:
    """
    Per transcript: (private start positions, total start positions).

    A start position is private when the resulting genomic footprint cannot be
    produced by any other transcript of the gene.  The ratio is the share of a
    molecule that actually carries evidence of its own identity.

    A transcript can be structurally identifiable (rank full) and still have a
    private fraction near zero: the distinguishing read exists but almost
    nothing lands on it.  That is the regime where quantifiers return confident
    answers with no support, and it is invisible to a rank test.
    """
    K = gene.n_iso
    priv = [0] * K
    tot = [0] * K
    for _sig, who in signature_owners(gene, R, F).items():
        solo = len(who) == 1
        for j, n in who.items():
            tot[j] += n
            if solo:
                priv[j] += n
    return list(zip(priv, tot))
