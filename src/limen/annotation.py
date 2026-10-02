"""
Annotation loading: GFF3 and GTF, plain or gzipped, standard library only.

Deliberately dependency free, so Level 0 runs anywhere Python does: no conda
requirement, no gffutils, no database build step.

Exon structure is taken from `exon` features when they exist.  Many
annotations (Phytozome among them) omit `exon` entirely and record only CDS
and UTR features, which together tile the mature transcript; those are used
as a fallback, per transcript, and abutting pieces are merged.  Which route
each transcript took is reported, because a parser that silently finds
nothing is worse than one that fails loudly.

Reference files are never modified.  Unwanted records (spike-ins such as
r_luc, organelles, unplaced contigs, transposable elements) are removed at
load time, so the same raw annotation stays usable by every other project.
"""
from __future__ import annotations

import gzip
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .identify import Gene, Transcript

__all__ = ["load_annotation", "detect_format", "AnnotationStats"]

_GTF_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')
_DEFAULT_EXCLUDE = ("r_luc",)
_FALLBACK_TYPES = ("CDS", "five_prime_UTR", "three_prime_UTR",
                   "5UTR", "3UTR", "UTR")


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def detect_format(path: str) -> str:
    """Return 'gtf' or 'gff3' by inspecting the first real attribute column."""
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attr = f[8]
            if "gene_id " in attr or "transcript_id " in attr:
                return "gtf"
            if "=" in attr:
                return "gff3"
    raise ValueError(f"cannot determine annotation format: {path}")


def _gff3_attrs(attr: str) -> Dict[str, str]:
    out = {}
    for part in attr.rstrip(";").split(";"):
        if "=" in part:
            k, _, v = part.partition("=")
            out[k.strip()] = v.strip()
    return out


def _gtf_attrs(attr: str) -> Dict[str, str]:
    return {k: v for k, v in _GTF_ATTR.findall(attr)}


class AnnotationStats:
    """What the parser actually found.  Print it; do not trust it silently."""

    def __init__(self) -> None:
        self.lines = 0
        self.exon_features = 0
        self.fallback_features = 0
        self.tx_from_exon = 0
        self.tx_from_fallback = 0
        self.skipped_seqid = 0
        self.skipped_excluded = 0
        self.skipped_te = 0
        self.orphan = 0
        self.genes = 0
        self.transcripts = 0
        self.seqids: Dict[str, int] = defaultdict(int)
        self.feature_types: Dict[str, int] = defaultdict(int)

    def __str__(self) -> str:
        top = sorted(self.seqids.items(), key=lambda kv: -kv[1])[:8]
        chrom = ", ".join(f"{k}({v:,})" for k, v in top)
        types = ", ".join(f"{k}={v:,}" for k, v in
                          sorted(self.feature_types.items(), key=lambda kv: -kv[1])[:8])
        route = (f"exon features {self.tx_from_exon:,} tx, "
                 f"CDS/UTR fallback {self.tx_from_fallback:,} tx")
        return (
            f"  lines read        : {self.lines:,}\n"
            f"  feature types     : {types}\n"
            f"  structure source  : {route}\n"
            f"  genes             : {self.genes:,}\n"
            f"  transcripts       : {self.transcripts:,}\n"
            f"  skipped (seqid)   : {self.skipped_seqid:,}\n"
            f"  skipped (exclude) : {self.skipped_excluded:,}\n"
            f"  skipped (TE)      : {self.skipped_te:,}\n"
            f"  orphan features   : {self.orphan:,}\n"
            f"  sequences         : {chrom}"
        )


def load_annotation(
    path: str,
    fmt: Optional[str] = None,
    seqids: Optional[Sequence[str]] = None,
    exclude: Sequence[str] = _DEFAULT_EXCLUDE,
    min_isoforms: int = 1,
    exclude_te: bool = True,
    verbose: bool = False,
) -> Tuple[List[Gene], AnnotationStats]:
    """
    Load genes from a GFF3 or GTF file.

    seqids       keep only these sequences (for example chromosome_01 ...)
    exclude      drop records whose seqid or identifiers contain these
                 substrings; defaults to the r_luc spike-in
    min_isoforms keep only genes with at least this many distinct transcripts
    exclude_te   drop transposable_element_gene entries
    """
    fmt = fmt or detect_format(path)
    keep = set(seqids) if seqids else None
    bad = tuple(e.lower() for e in exclude)
    st = AnnotationStats()

    exon_parts: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    fallback_parts: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    tx_meta: Dict[str, Tuple[str, str]] = {}
    child_parent: Dict[str, str] = {}
    te_ids: Set[str] = set()

    with _open(path) as fh:
        for line in fh:
            st.lines += 1
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqid, ftype, start, end, strand = f[0], f[2], f[3], f[4], f[6]
            attr = f[8]
            st.feature_types[ftype] += 1

            if fmt == "gff3":
                a = _gff3_attrs(attr)
                if "ID" in a:
                    if "Parent" in a:
                        child_parent[a["ID"]] = a["Parent"].split(",")[0]
                    if ftype == "transposable_element_gene":
                        te_ids.add(a["ID"])
                tid = a.get("Parent", "").split(",")[0] or a.get("transcript_id", "")
            else:
                a = _gtf_attrs(attr)
                tid = a.get("transcript_id", "")

            is_exon = ftype == "exon"
            is_fallback = ftype in _FALLBACK_TYPES
            if not (is_exon or is_fallback):
                continue
            if not tid:
                st.orphan += 1
                continue

            low = (seqid + " " + tid + " " + attr).lower()
            if any(b in low for b in bad):
                st.skipped_excluded += 1
                continue
            if keep is not None and seqid not in keep:
                st.skipped_seqid += 1
                continue

            # GFF and GTF are 1-based inclusive; convert to 0-based half-open
            span = (int(start) - 1, int(end))
            if is_exon:
                exon_parts[tid].append(span)
                st.exon_features += 1
            else:
                fallback_parts[tid].append(span)
                st.fallback_features += 1
            tx_meta.setdefault(tid, (seqid, strand))
            st.seqids[seqid] += 1

            if fmt == "gtf":
                child_parent.setdefault(tid, a.get("gene_id", tid))

    def _gene_of(tid: str) -> str:
        gid, hops = tid, 0
        while gid in child_parent and hops < 10:
            gid = child_parent[gid]
            hops += 1
        return gid

    def _merge(parts: List[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
        parts.sort()
        out: List[Tuple[int, int]] = []
        for s, e in parts:
            if out and s <= out[-1][1]:        # abutting or overlapping: one exon
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        return tuple(out)

    by_gene: Dict[str, List[Transcript]] = defaultdict(list)
    # sorted, not set order: hash randomisation would otherwise make
    # which transcript survives deduplication differ between runs
    for tid in sorted(set(exon_parts) | set(fallback_parts)):
        if tid in exon_parts:                  # prefer real exon features
            structure = _merge(exon_parts[tid])
            st.tx_from_exon += 1
        else:
            structure = _merge(fallback_parts[tid])
            st.tx_from_fallback += 1
        gid = _gene_of(tid)
        if exclude_te and gid in te_ids:
            st.skipped_te += 1
            continue
        by_gene[gid].append(Transcript(tid, structure))

    genes: List[Gene] = []
    for gid, txs in by_gene.items():
        uniq: Dict[Tuple, Transcript] = {}
        for t in txs:                          # identical structures are one unit
            uniq.setdefault(t.exons, t)
        if len(uniq) < min_isoforms:
            continue
        seqid, strand = tx_meta.get(txs[0].id, ("", "."))
        genes.append(Gene(gid, tuple(sorted(uniq.values(), key=lambda t: t.id)),
                          seqid, strand))

    genes.sort(key=lambda g: g.id)
    st.genes = len(genes)
    st.transcripts = sum(g.n_iso for g in genes)
    if verbose:
        print(f"[limen] {path} ({fmt})\n{st}", file=sys.stderr)
    return genes, st
