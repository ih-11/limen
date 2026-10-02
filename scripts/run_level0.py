#!/usr/bin/env python3
"""
Level 0 sweep for one annotated transcriptome.

Produces, for a ladder of read lengths:

  <prefix>.genes.tsv        per gene: rank of the incidence matrix, whether
                            the gene is structurally identifiable, and by how
                            many dimensions it falls short if not
  <prefix>.transcripts.tsv  per transcript: how many of its possible reads are
                            diagnostic (no sibling transcript can produce the
                            same genomic footprint), and the same quantity
                            expressed per megabase of sequencing
  <prefix>.run.json         exactly which file and parameters produced the
                            above, so a result can always be traced back

Reads no sequence.  An annotation is the only input.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import limen
from limen import diagnostic_fraction, identifiability, load_annotation

DEFAULT_LADDER = "75,100,150,250,300,500,1000,2000,4000,8000,16000"


def _sweep_gene(args):
    """One gene across the whole ladder. Top level so it can be pickled."""
    gene, reads, frags = args
    grows, trows = [], []
    for R in reads:
        for F in frags:
            if F is not None and F < R:
                continue
            r = identifiability(gene, R, F)
            grows.append((gene.id, gene.seqid, gene.n_iso, R,
                          "NA" if F is None else F, r.rank, r.deficiency,
                          int(r.identifiable), r.n_patterns))
            for t, (diag, tot) in zip(gene.transcripts,
                                      diagnostic_fraction(gene, R, F)):
                frac = diag / tot if tot else 0.0
                # a read costs min(read length, transcript length) bases
                cost = min(R, t.length) or 1
                trows.append((gene.id, t.id, t.length, gene.n_iso, R,
                              "NA" if F is None else F, diag, tot,
                              f"{frac:.6f}", f"{frac * 1e6 / cost:.3f}"))
    return grows, trows


def _digest(path: Path, cap: int = 64 << 20) -> str:
    """md5 of up to the first 64 MB, enough to catch a changed annotation."""
    h, n = hashlib.md5(), 0
    with open(path, "rb") as fh:
        while n < cap:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
            n += len(chunk)
    return f"md5:{h.hexdigest()}" + ("" if n < cap else f" (first {cap} bytes)")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Structural identifiability and diagnostic evidence "
                    "across read lengths, from an annotation alone.")
    p.add_argument("annotation", help="GFF3 or GTF, optionally gzipped")
    p.add_argument("--species", required=True,
                   help="label used for output filenames and the sidecar")
    p.add_argument("--outdir", default=None,
                   help="default: $LIMEN_WORK/results/level0")
    p.add_argument("--read-lengths", default=DEFAULT_LADDER)
    p.add_argument("--fragment-lengths", default=None,
                   help="comma separated; enables paired-end mode")
    p.add_argument("--min-isoforms", type=int, default=2,
                   help="single-isoform genes are trivially identifiable "
                        "[%(default)s]")
    p.add_argument("--max-isoforms", type=int, default=0,
                   help="skip genes above this many isoforms, 0 = no limit; "
                        "cost grows quickly with isoform count")
    p.add_argument("--seqids", default=None,
                   help="comma separated; keep only these sequences")
    p.add_argument("--exclude", default="r_luc",
                   help="comma separated substrings to drop [%(default)s]")
    p.add_argument("--keep-te", action="store_true",
                   help="keep transposable_element_gene entries")
    p.add_argument("--limit", type=int, default=0,
                   help="use only the first N genes, for a quick trial")
    p.add_argument("--jobs", type=int, default=1,
                   help="parallel worker processes [%(default)s]")
    arg = p.parse_args(argv)

    reads = [int(x) for x in arg.read_lengths.split(",") if x]
    frags = ([int(x) for x in arg.fragment_lengths.split(",")]
             if arg.fragment_lengths else [None])
    outdir = Path(arg.outdir or
                  Path(os.environ.get("LIMEN_WORK", Path.home() / "work/limen"))
                  / "results" / "level0")
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / arg.species
    ann = Path(arg.annotation)

    t_start = time.time()
    genes, stats = load_annotation(
        str(ann),
        seqids=arg.seqids.split(",") if arg.seqids else None,
        exclude=[e for e in arg.exclude.split(",") if e],
        exclude_te=not arg.keep_te,
        verbose=True)
    t_load = time.time() - t_start

    kept = [g for g in genes
            if g.n_iso >= arg.min_isoforms
            and (arg.max_isoforms == 0 or g.n_iso <= arg.max_isoforms)]
    skipped_big = sum(1 for g in genes
                      if arg.max_isoforms and g.n_iso > arg.max_isoforms)
    if arg.limit:
        kept = kept[:arg.limit]

    print(f"[limen] {len(kept):,} genes enter the sweep "
          f"({len(genes):,} loaded, {skipped_big:,} over the isoform cap)",
          file=sys.stderr)
    if not kept:
        print("[limen] no genes meet the isoform threshold; writing empty "
              "tables and the sidecar, which is itself the result",
              file=sys.stderr)

    work = [(g, reads, frags) for g in kept]
    t0 = time.time()
    results = []
    if arg.jobs > 1:
        import multiprocessing as mp
        with mp.Pool(arg.jobs) as pool:
            for i, out in enumerate(pool.imap_unordered(_sweep_gene, work, 64), 1):
                results.append(out)
                if i % 2000 == 0 or i == len(work):
                    el = time.time() - t0
                    print(f"[limen]   {i:,}/{len(work):,} genes  {el:,.0f}s  "
                          f"eta {el/i*(len(work)-i):,.0f}s", file=sys.stderr)
    else:
        for i, w in enumerate(work, 1):
            results.append(_sweep_gene(w))
            if i % 2000 == 0 or i == len(work):
                el = time.time() - t0
                print(f"[limen]   {i:,}/{len(work):,} genes  {el:,.0f}s  "
                      f"eta {el/i*(len(work)-i):,.0f}s", file=sys.stderr)
    t_sweep = time.time() - t0

    gpath, tpath = f"{prefix}.genes.tsv", f"{prefix}.transcripts.tsv"
    ng = nt = 0
    with open(gpath, "w", newline="") as gf, open(tpath, "w", newline="") as tf:
        gw = csv.writer(gf, delimiter="\t", lineterminator="\n")
        tw = csv.writer(tf, delimiter="\t", lineterminator="\n")
        gw.writerow(["gene_id", "seqid", "n_isoforms", "read_length",
                     "fragment_length", "rank", "deficiency", "identifiable",
                     "n_patterns"])
        tw.writerow(["gene_id", "tx_id", "tx_len", "n_isoforms", "read_length",
                     "fragment_length", "diagnostic", "total", "frac",
                     "diag_per_mb"])
        for grows, trows in results:
            gw.writerows(grows); ng += len(grows)
            tw.writerows(trows); nt += len(trows)

    sidecar = {
        "species": arg.species,
        "annotation": {
            "path": str(ann.resolve()),
            "name": ann.name,
            "bytes": ann.stat().st_size,
            "modified": datetime.fromtimestamp(ann.stat().st_mtime,
                                               timezone.utc).isoformat(),
            "digest": _digest(ann),
        },
        "parameters": {
            "read_lengths": reads,
            "fragment_lengths": frags,
            "min_isoforms": arg.min_isoforms,
            "max_isoforms": arg.max_isoforms,
            "seqids": arg.seqids,
            "exclude": arg.exclude,
            "exclude_te": not arg.keep_te,
            "limit": arg.limit,
        },
        "parse": {
            "genes_loaded": stats.genes,
            "transcripts_loaded": stats.transcripts,
            "genes_swept": len(kept),
            "skipped_over_isoform_cap": skipped_big,
            "structure_from_exon_features": stats.tx_from_exon,
            "structure_from_cds_utr_fallback": stats.tx_from_fallback,
            "skipped_seqid": stats.skipped_seqid,
            "skipped_excluded": stats.skipped_excluded,
            "skipped_te": stats.skipped_te,
        },
        "output": {"genes_tsv": gpath, "gene_rows": ng,
                   "transcripts_tsv": tpath, "transcript_rows": nt},
        "timing_seconds": {"load": round(t_load, 1), "sweep": round(t_sweep, 1)},
        "limen_version": limen.__version__,
        "python": sys.version.split()[0],
        "run_at": datetime.now(timezone.utc).isoformat(),
    }
    with open(f"{prefix}.run.json", "w") as fh:
        json.dump(sidecar, fh, indent=2)

    print(f"[limen] wrote {ng:,} gene rows and {nt:,} transcript rows "
          f"in {t_sweep:,.0f}s", file=sys.stderr)
    print(f"[limen] {prefix}.{{genes,transcripts}}.tsv  +  .run.json",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
