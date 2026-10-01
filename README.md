# limen

**Where a molecular distinction becomes visible to sequencing.**

*limen* (Latin, *threshold*; in psychophysics, the smallest stimulus that can
still be detected) measures the read length at which two different RNA
molecules stop being confusable — computed from an annotation alone, before
any sequencing is done.

## The question

The molecule in the cell is not the molecule in the tube, is not the read, is
not the FASTQ, is not the count table. Every step loses something, and in a
real experiment the starting state is never observed, so the loss cannot be
measured. This project defines the starting state and measures the loss.

Concretely: **what can long reads see that short reads cannot?** Stated
mathematically that is a question of *identifiability* — when do two genuinely
different transcript populations produce byte-identical data?

### Worked example

A gene `A–[B]–M(2000 nt)–[C]–D` where B and C are independently skippable:
four isoforms, two splicing decisions ~2000 nt apart.

| population | composition |
|---|---|
| 1 | 50% keep both + 50% skip both |
| 2 | 50% keep B only + 50% keep C only |

Different molecules. At 150 bp the observed data is **identical** — not noisy,
not underpowered: the same numbers. Short reads see two marginals and four
unknowns. At 2002 nt — the shortest read reaching from inside B to inside C —
they separate.

Quantifiers still return four abundances. In that regime the EM picks one
point from an infinite family of equally valid solutions.

## Method

Every possible read contributes one row of a 0/1 incidence matrix **M**: a 1
for each transcript that could have produced it.

- `rank(M) == K` → abundances uniquely determined. **Identifiable.**
- `rank(M) < K` → different mixtures give identical data. `K - rank(M)` is the
  number of invisible dimensions.

This is the infinite-depth, error-free **ceiling**. Real pipelines land below
it.

## Levels

| | measures | needs |
|---|---|---|
| **0** identifiability | the ceiling vs read length | annotation only |
| **1** simulation | what real tools recover; the gap | + genome FASTA, aligner, quantifier |
| **2** real data | whether the predicted gap holds | matched long/short-read libraries |

Level 0 never reads sequence. A GFF3 or GTF is the only input.

## Setup

Requires Python ≥ 3.9 and numpy. Nothing else.

```bash
git clone https://github.com/ih-11/limen.git
cd limen
pip install -e ".[dev]"
pytest -q
```

Activate whichever environment you use first — the repo does not care which:

```bash
conda activate ih      # macOS
conda activate ibnu    # WSL
```

### Paths

No absolute path appears anywhere in this repository. Point it at your data
with two environment variables, set in `~/.zshrc` (macOS) or `~/.bashrc`
(WSL):

```bash
# macOS
export LIMEN_REF="$HOME/Code/ReferenceGenome"
export LIMEN_WORK="$HOME/work/limen"

# WSL
export LIMEN_REF="/mnt/f/RA/Downstream/Project4_LIMEN/ReferenceGenome"
export LIMEN_WORK="$HOME/work/limen"
```

Code lives in git, large data never does, and generated output goes to
`$LIMEN_WORK`. On WSL keep working data under `$HOME` rather than `/mnt/`,
which is noticeably slower.

## Prior work

The identifiability framework for short-read RNA-seq is due to Hiller *et al.*
(2009) *Bioinformatics* **25**:3056, and Ferrer-Bonsoms *et al.* (2022)
*Bioinformatics* **38**:1491, the latter computing identifiability as a
function of read and fragment length for human only, 75–300 bp, assuming
infinite depth, and not treating long reads.

This is an independent Python implementation written from the published
descriptions; no code is derived from those repositories. It extends the
question to the long-read regime and across transcriptomes of differing
architecture.

## Licence

MIT.
