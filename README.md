# limen

**Measuring what sequencing can know, and what it only appears to know.**

*limen* (Latin, *threshold*; in psychophysics, the smallest stimulus that can
still be detected).

## The problem

Sequencing does not observe the transcriptome. It observes whatever survives
extraction, fragmentation, reverse transcription, amplification, basecalling,
alignment and statistical inference. The molecule in the cell is not the
molecule in the tube, is not the read, is not the FASTQ, is not the count
table.

Each of those steps discards information. None of them can be measured in a
real experiment, because the starting state is never observed. We compare the
end of the chain against the end of another chain and call the difference a
result.

This project defines the starting state, so the loss becomes measurable.

## The question

Long-read sequencing is usually chosen by category: *I study structural
variants, therefore long reads.* The choice is rarely quantitative. The
question here is the quantitative one:

> What can long reads know that short reads cannot, and to what extent?

Not "which platform is better", which depends on the question asked, but where
exactly the boundary sits, what sets its position, and when the short-read
answer stops being merely blurry and becomes actively wrong.

## Approach

Three levels, each a different claim.

| | asks | needs |
|---|---|---|
| **0** | what is knowable in principle, at a given read length | annotation only |
| **1** | what a real pipeline actually recovers | simulated reads, aligner, quantifier |
| **2** | whether the predicted behaviour holds on real libraries | matched long and short read data |

Level 0 gives a ceiling: infinite depth, no errors, perfect annotation. Real
pipelines land below it.

**The distance between the ceiling and what is achieved is the result.** Level
0 alone is a known quantity in one organism. The gap is not.

Two further axes that no existing study covers:

- **Transcriptome architecture as a variable.** The relationship between read
  length and recoverable information depends on how exons are sized and
  spaced, and every published analysis is human. Six genomes spanning
  intron-poor yeasts, intron-dense green algae, and an allotetraploid plant
  make architecture an independent variable rather than a constant.
- **Reads generated under full control.** Fragment length, error profile,
  truncation, chimera rate, depth and annotation completeness become separate
  dials rather than bundled platform presets, so their contributions can be
  told apart instead of attributed wholesale to "read length".

## A worked example

A gene `A-[B]-M(2000 nt)-[C]-D`, where B and C are independently skippable:
four isoforms, two splicing decisions about 2000 nt apart in the mature
message.

| population | composition |
|---|---|
| 1 | 50% keep both, 50% skip both |
| 2 | 50% keep B only, 50% keep C only |

Different molecules, possibly different proteins. At 150 bp the observed data
from these two populations is *identical*. Not noisy, not underpowered: the
same numbers. Short reads see two marginals (how often B appears, how often C
appears) against four unknowns.

At 2002 nt, the shortest read reaching from inside B to inside C, they
separate cleanly.

Quantifiers are asked for four abundances and return four regardless. In the
deficient regime the EM selects one point from an infinite family of equally
valid solutions, guided by its priors rather than by the data. Those invented
answers are detectable: Chen *et al.* (2025) found short-read-specific "major
isoforms", tested thirteen by digital PCR, and showed they were not the
dominant molecules.

## Level 0: how the ceiling is computed

Each possible read contributes one row of a 0/1 incidence matrix **M**, with a
1 for every transcript that could have produced it.

- `rank(M) == K`: abundances are uniquely determined. Identifiable.
- `rank(M) < K`: different transcript mixtures produce identical data, and no
  algorithm, no sequencing depth and no replication can separate them.
  `K - rank(M)` counts the invisible dimensions.

Identifiability is a standard tool, not a contribution of this project (see
Prior work). It is used here because it is the correct way to state the
ceiling. Level 0 never reads sequence: a GFF3 or GTF is the only input.

## Status

Early. Level 0 is being built and validated. Levels 1 and 2 are not started.

## Setup

Python 3.9 or newer, and numpy. Nothing else.

```bash
git clone https://github.com/ih-11/limen.git
cd limen
pip install -e ".[dev]"
pytest -q
```

Activate whichever environment you use first. The repository does not care
which, and never names one:

```bash
conda activate ih      # macOS
conda activate ibnu    # WSL
```

### Paths

No absolute path appears in any committed file. Point the code at your data
with two environment variables, in `~/.zshrc` on macOS or `~/.bashrc` on WSL:

```bash
# macOS
export LIMEN_REF="$HOME/Code/ReferenceGenome"
export LIMEN_WORK="$HOME/work/limen"

# WSL
export LIMEN_REF="/mnt/f/RA/Downstream/Project4_LIMEN/ReferenceGenome"
export LIMEN_WORK="$HOME/work/limen"
```

Code lives in git. Large data never does. Generated output goes to
`$LIMEN_WORK`. On WSL, keep working data under `$HOME` rather than under
`/mnt/`, which is markedly slower.

Reference files are never modified. Spike-ins, organelles and unplaced
contigs are filtered when the annotation is read.

## Prior work

Structural identifiability of isoform deconvolution was established for
short-read RNA-seq by:

- Hiller *et al.* (2009) *Bioinformatics* **25**:3056. Identifiability
  conditions for isoform deconvolution. Read length is not a variable.
- Ferrer-Bonsoms *et al.* (2022) *Bioinformatics* **38**:1491. Identifiability
  as a function of read and fragment length, used to select a library. Human
  only, 75 to 300 bp, infinite depth assumed, long reads not considered.

This repository contains an independent Python implementation written from the
published descriptions. No code is derived from those projects. Reproducing
their human results is used as a correctness check on the implementation, not
as a result.

What is new here is the extension of the read-length axis into the long-read
range, the treatment of transcriptome architecture as an independent variable,
and the measurement of the distance between the structural ceiling and what
real pipelines recover.

## Licence

MIT.
