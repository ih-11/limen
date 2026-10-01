# limen

**Measuring what sequencing can know, and what it only appears to know.**

*limen* (Latin, *threshold*; in psychophysics, the smallest stimulus that can
still be detected). Given an annotation, it computes the read length and
sequencing depth at which two different RNA molecules stop being confusable,
before any sequencing is done.

---

## The problem

Sequencing does not observe the transcriptome. It observes whatever survives
extraction, fragmentation, reverse transcription, amplification, basecalling,
alignment and statistical inference. The molecule in the cell is not the
molecule in the tube, is not the read, is not the FASTQ, is not the count
table.

Each step discards information, and none of the losses can be measured in a
real experiment, because the starting state is never observed. We compare the
end of one chain against the end of another and call the difference a result.

This project defines the starting state, so the loss becomes measurable.

## The question

Long-read sequencing is usually chosen by category: *I study structural
variants, therefore long reads*. The choice is rarely quantitative. The
question here is the quantitative one:

> What can long reads know that short reads cannot, and to what extent?

Not which platform is better, which depends on the question asked, but where
the boundary sits, what sets its position, and when a short-read answer stops
being merely imprecise and becomes unsupported.

## Approach

| | asks | needs | status |
|---|---|---|---|
| **Level 0** | what is knowable in principle at a given read length | annotation only | implemented |
| **Level 1** | what a real pipeline actually recovers | simulated reads, aligner, quantifier | not started |
| **Level 2** | whether the predicted behaviour holds on real libraries | matched long and short read data | not started |

**The distance between the ceiling and what is achieved is the result.**
Level 0 alone is a known quantity in one organism. The gap is not.

Two axes no existing study covers:

- **Transcriptome architecture as a variable.** How read length relates to
  recoverable information depends on how exons are sized and spaced, and every
  published analysis is human. A panel spanning intron-poor fungi, intron-rich
  plants and algae, and a polyploid genome makes architecture an independent
  variable rather than a constant.
- **Reads generated under full control.** Fragment length, error profile,
  truncation, chimera rate, depth and annotation completeness become separate
  dials rather than bundled platform presets, so their contributions can be
  told apart instead of attributed wholesale to read length.

## What Level 0 computes

Two quantities, which answer different questions.

**Identifiability (the ceiling).** Every possible read contributes one row of a
0/1 incidence matrix `M`, with a 1 for each transcript that could have produced
it. If `rank(M) == K` the transcript abundances are uniquely determined. If
`rank(M) < K`, genuinely different mixtures produce identical data and no
algorithm, depth or replication can separate them. `K - rank(M)` counts the
invisible dimensions.

This assumes unlimited depth, so it reports whether a distinguishing read
*exists*, not whether it would ever be seen.

**Diagnostic evidence (the conditioning).** For each transcript, the share of
its possible reads whose genomic footprint no sibling transcript could produce.
Such a read is direct evidence that the molecule was present; a read a sibling
could also explain is real data but silent about its origin.

Zero does not mean a transcript gets no reads. It means every read it can
produce is also explainable by a sibling, so its abundance is reachable only by
subtracting the others. That is qualitatively different from a small fraction:
small needs more depth, zero cannot be fixed by depth at all.

A transcript can be structurally identifiable and diagnostic nowhere. A rank
test cannot see that, and it is the regime in which quantifiers return
confident estimates with no support.

Level 0 never reads sequence. A GFF3 or GTF is the only input.

## Status

Level 0 is implemented and validated. Levels 1 and 2 are not started.

Implemented: identifiability by matrix rank, diagnostic fraction per
transcript, single-end and paired-end, GFF3 and GTF reading with exon
reconstruction from CDS and UTR features where `exon` is absent.

Pending: reproduction of published human short-read identifiability rates as an
external check on the implementation, the architecture sweep across the full
genome panel, and all of Levels 1 and 2.

## Install

Python 3.9 or newer, and numpy. Nothing else.

```bash
git clone https://github.com/ih-11/limen.git
cd limen
pip install -e ".[dev]"
pytest -q
```

Activate whichever environment you use first. The repository does not care
which, and never names one.

### Paths

No absolute path appears in any committed file. Point the code at your data
with two environment variables, in `~/.zshrc` on macOS or `~/.bashrc` on WSL:

```bash
export LIMEN_REF="$HOME/Code/ReferenceGenome"      # annotations and genomes
export LIMEN_WORK="$HOME/work/limen"               # generated output
```

Code lives in git. Large data never does. Generated output goes to
`$LIMEN_WORK`. On WSL, keep working data under `$HOME` rather than under
`/mnt/`, which is markedly slower.

Reference files are never modified. Spike-ins, organelles, unplaced contigs
and transposable element genes are filtered when the annotation is read.

## Use

```python
from limen import load_annotation, identifiability, diagnostic_fraction

genes, stats = load_annotation("annotation.gff3", verbose=True)

gene = genes[0]
r = identifiability(gene, 150)
print(r.rank, "of", r.n_iso, "identifiable" if r.identifiable else "deficient")

for tx, (diag, total) in zip(gene.transcripts, diagnostic_fraction(gene, 150)):
    print(tx.id, f"{100 * diag / total:.1f}% diagnostic")
```

When comparing read lengths, express results per unit of sequencing effort
rather than per read. A read costs `min(read length, transcript length)` bases,
and comparing per read overstates the long-read advantage by roughly the ratio
of the read lengths.

## Validation

Correctness is established three ways, in increasing order of strength.

1. **Toy genes with hand-derived answers.** The central case is four isoforms
   produced by two independent exon skips about 2000 nt apart: rank 3 of 4
   below 2002 nt, rank 4 at and above it, where 2002 is the shortest read
   reaching from inside one skipped exon to inside the other.
2. **Brute-force cross-checks on randomised genes.** The fast segmentation
   algorithms are compared against naive enumeration over thousands of random
   gene structures, single-end and paired-end, on every test run. This is not
   decorative: it caught three separate bugs during development, each of which
   produced plausible but wrong numbers that the toy genes alone did not expose.
3. **Published results.** Reproducing the human short-read identifiability
   rates reported by Ferrer-Bonsoms et al. (2022) is treated as a prerequisite
   before the read-length axis is extended beyond their range. It is a
   correctness check on the implementation, not a result.

## Layout

```
src/limen/   identify.py     the matrix, its rank, and diagnostic evidence
             annotation.py   GFF3 and GTF, standard library only
             __init__.py
test/        toy expectations and brute-force cross-checks
scripts/     batch jobs: produce data, write TSV, no plots
notebook/    interpret data, make figures, decide. reads what scripts wrote
ref/         manifest of genome sources. never the genomes themselves
docs/        METHODS.md, RESEARCH_PLAN.md
```

The split matters: if it must be correct it goes in `src/` with a test, if it
must be looked at it goes in a notebook, and anything heavy enough to run
unattended goes in `scripts/` and writes a file the notebook then reads.

## Prior work

Structural identifiability of isoform deconvolution was established for
short-read RNA-seq by Hiller et al. (2009), who derived the conditions but did
not treat read length as a variable, and by Ferrer-Bonsoms et al. (2022), who
computed identifiability as a function of read and fragment length in order to
select a library. The latter covers the human transcriptome only, read lengths
of 75 to 300 nt, assumes unlimited depth, and does not consider long reads.

This repository contains an independent Python implementation written from the
published descriptions. No code is derived from those projects.

What is new here is the extension of the read-length axis into the long-read
range, the treatment of transcriptome architecture as an independent variable,
the replacement of binary identifiability with a measure of how much evidence
exists, and the measurement of the distance between the structural ceiling and
what real pipelines recover.

The research plan, with full citations, is in `docs/RESEARCH_PLAN.md`.

## Licence

MIT.
