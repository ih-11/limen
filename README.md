# limen

**Measuring what sequencing can know, and what it only appears to know.**

*limen* (Latin, *threshold*; in psychophysics, the smallest stimulus that can
still be detected). Given an annotation, it computes the read length and the
sequencing depth at which two different RNA molecules stop being confusable,
before any sequencing has been done.

This README is written to be self-contained. Anyone picking the project up,
including a future conversation with no memory of the earlier ones, should be
able to read this file and know what the project is, what has been built, what
has been found, what is already known to be wrong, and what comes next.

---

## 1. The problem

Sequencing does not observe the transcriptome. It observes whatever survives
extraction, fragmentation, reverse transcription, amplification, basecalling,
alignment and statistical inference. The molecule in the cell is not the
molecule in the tube, is not the read, is not the FASTQ, is not the count
table.

Every one of those steps discards information, and in a real experiment none
of the losses can be measured, because the starting state is never observed.
What the field does instead is compare the end of one transformation chain
against the end of another and attribute the difference to the technologies
rather than to the chains.

This project defines the starting state, so that the loss becomes measurable.

## 2. The question

Long-read sequencing is usually chosen by category of application. A
researcher studying structural variants uses long reads because that is what
one does for structural variants. The choice is rarely quantitative, and the
justification is usually a benchmark performed in a different organism with a
different transcriptome.

The question here is the quantitative one:

> What can long reads know that short reads cannot, and to what extent?

Not which platform is better, which depends entirely on the question being
asked, but where the boundary sits, what determines its position, and when a
short-read answer stops being merely imprecise and becomes unsupported by the
data.

## 3. Approach: three levels

| | asks | needs | status |
|---|---|---|---|
| **Level 0** | what is knowable in principle at a given read length | annotation only | implemented and running |
| **Level 1** | what a real pipeline actually recovers | simulated reads, aligner, quantifier | not started |
| **Level 2** | whether the predicted behaviour holds on real libraries | matched long and short read data | not started |

Level 0 computes a ceiling by assuming unlimited depth, no sequencing error
and a complete annotation. Real pipelines land below that ceiling.

**The distance between the ceiling and what is actually achieved is the
result.** Level 0 on its own is a known quantity in one organism, already
published for human short reads. The gap is not known for anything.

## 4. What Level 0 computes

Two quantities, which answer different questions and disagree in an
informative way.

### Identifiability, which is the ceiling

Every possible read contributes one row of a 0/1 incidence matrix `M`, with a
1 in the column of each transcript that could have produced that read. Entries
are 0 or 1 rather than counts because a genomic footprint fixes the start
position, which fixes the position within any transcript able to produce it.

- `rank(M) == K`: the transcript abundances are uniquely determined by the
  data. The gene is identifiable.
- `rank(M) < K`: genuinely different mixtures of transcripts produce
  byte-identical data, and no algorithm, no sequencing depth and no number of
  replicates can separate them. `K - rank(M)` counts the invisible dimensions.

### Diagnostic evidence, which is the conditioning

Rank asks whether a distinguishing read exists somewhere. It does not ask
whether that read would ever be observed. A gene counts as identifiable even
if the only thing separating two of its transcripts is a twenty nucleotide
window that almost no read will ever cover.

So for each transcript we also compute the share of its possible reads whose
genomic footprint no sibling transcript could produce. Such a read is direct
evidence that this particular molecule was present. A read that a sibling
could also explain is real data, but it is silent about which molecule it came
from.

A diagnostic fraction of zero does not mean a transcript receives no reads. It
means that every read it can produce is also explainable by a sibling, so its
abundance is reachable only by subtracting the others. That is algebraically
valid and statistically unstable, and it is qualitatively different from a
small fraction. A small fraction needs more depth. Zero cannot be fixed by
depth at all, because zero multiplied by any number of reads is still zero.

A transcript can therefore be structurally identifiable and diagnostic
nowhere. A rank test cannot see this, and it is precisely the regime in which
quantification methods return confident estimates with no supporting evidence.

Level 0 never reads sequence. A GFF3 or GTF is the only input.

## 5. Findings so far

These are preliminary, from Level 0 only, and the scope of each is stated
because the numbers are meaningless without it.

### The ceiling is high

On the 300 most isoform-rich genes of *Chlamydomonas reinhardtii*, meaning the
hardest genes in that genome, **97.7 per cent are structurally identifiable at
a 100 nucleotide read length**, rising to 100 per cent by 5 kb. Long reads buy
roughly two percentage points.

The mechanism is that annotated isoforms almost always differ by at least one
local feature. In a gene with twenty-one isoforms and twenty-three exons, each
isoform tends to carry at least one junction or exon boundary unique to it,
which produces a matrix row with a single 1 in it, and enough such rows make a
full-rank submatrix. The pathological case, in which isoforms differ only in
how distant decisions are combined, turns out to be rare.

Whatever long reads buy, it is therefore not resolving-in-principle given a
complete annotation.

### But a quarter of transcripts carry no evidence of their own identity

In the same gene set, at a 150 nucleotide read length:

| | |
|---|---|
| identifiable by rank | 98.0 per cent |
| transcripts with **zero** diagnostic reads | **25.3 per cent** |
| under 1 per cent diagnostic | 41.3 per cent |
| under 5 per cent diagnostic | 74.2 per cent |
| median diagnostic fraction | 1.7 per cent |

A quarter of transcripts in the most isoform-rich genes appear in a
quantifier's output only as the residue of subtracting their siblings.

This is the quantitative form of a phenomenon that has already been confirmed
experimentally. Chen and colleagues found transcripts reported as dominant by
short-read analysis, tested thirteen of them by digital PCR, and showed they
were not the dominant molecules in the sample.

### Per base sequenced, the advantage is about threefold and arrives as a step

Comparisons must be made per unit of sequencing effort rather than per read,
because a 4 kb read costs roughly twenty-six times a 150 nucleotide read. A
read also costs `min(read length, transcript length)` bases, never more, since
you cannot spend 8 kb of sequencing to read a 2.7 kb molecule.

On that basis, diagnostic reads per megabase run at roughly 115 across the
short-read range, **dip to about 89 at 1 kb**, and then climb to a plateau of
about 275 once reads clear whole molecules.

The dip is the interesting part. Between 150 nucleotides and 1 kb the
diagnostic fraction rises about fivefold while the cost per read rises nearly
sevenfold, so per base the result goes slightly backwards. Partial long reads
are the worst of both worlds: you pay for length and do not get the thing that
length is for. The benefit appears as a step function when reads clear the
transcript, not as a gradient.

### The architecture panel is confounded with annotation depth

Level 0 was run across five annotations. The isoform content is:

| species | genes | transcripts | isoforms per gene | multi-isoform genes |
|---|---|---|---|---|
| *Komagataella phaffii* | 5,040 | 5,040 | 1.00 | 0 |
| *Nicotiana benthamiana* | 59,814 | 59,814 | 1.00 | 0 |
| *Oryza sativa* | 37,858 | 44,714 | 1.18 | 5,384 |
| *Arabidopsis thaliana* | 33,016 | 54,369 | 1.65 | 11,058 |
| *Chlamydomonas reinhardtii* | 16,883 | 31,858 | 1.89 | 7,287 |

Two of the five annotations call no alternative isoforms whatsoever. The rice
file is the IRGSP *representative* transcript set, deliberately reduced to
approximately one model per gene, so its 1.18 is a property of which file is
on disk rather than of rice.

The consequence is serious and must not be glossed over. **The two species
with rich isoform annotation are also the two whose annotation projects set
out to call isoforms.** Transcriptome architecture and annotation effort are
completely confounded in this panel, so a claim of the form "species A needs
longer reads than species B" cannot be defended.

### The reframe this forces

Species is the wrong unit of analysis. Architecture should be measured per
gene, not per species.

Across *Chlamydomonas*, *Arabidopsis* and *Oryza* there are already **23,729
multi-isoform genes** spanning a wide range of exon counts, exon lengths,
transcript lengths and isoform multiplicities. The relationship between read
length and recoverable information can be measured directly across that range,
using gene-level architecture as the predictor. Species then becomes a
robustness check rather than the explanatory variable, and the confound stops
mattering because no claim rests on it.

This is a better design in any case. "The required read length is predictable
from gene structure" is both stronger and more useful than a ranking of
species, and it is the question a laboratory would actually ask.

## 6. Data

### Where things are

| | path |
|---|---|
| repository | `~/Code/limen` on both machines, `github.com/ih-11/limen` |
| annotations and genomes, WSL | `/mnt/f/RA/Downstream/Project4_LIMEN/ReferenceGenome` |
| annotations and genomes, Mac | `~/Code/ReferenceGenome` |
| generated output, both | `~/work/limen` |
| project archive, drafts, final figures | `/mnt/f/RA/Downstream/Project4_LIMEN` |

No absolute path appears in any committed file. The code is pointed at data
by two environment variables, set in `~/.bashrc` on WSL and `~/.zshrc` on the
Mac:

```bash
export LIMEN_REF="/mnt/f/RA/Downstream/Project4_LIMEN/ReferenceGenome"   # WSL
export LIMEN_REF="$HOME/Code/ReferenceGenome"                            # Mac
export LIMEN_WORK="$HOME/work/limen"                                     # both
```

Code lives in git. Large data never does. Generated output goes to
`$LIMEN_WORK` and is excluded by `.gitignore`.

### Annotations currently available

All carry a `std_r_luc` suffix, meaning a Renilla luciferase spike-in cassette
was added to the reference for a separate project. It is filtered out at load
time and the reference files themselves are never modified.

| species | file | source | note |
|---|---|---|---|
| *C. reinhardtii* | `CreinhardtiiCC_4532_707_v6.0.chr.std_r_luc.gff3` | Phytozome | no `exon` features, structure rebuilt from CDS and UTR |
| *A. thaliana* | `Araport11_GFF3_genes_transposons.201606.std_r_luc.gff3` | Araport | mixed, mostly `exon`, includes organelles `ChrM` and `ChrC` |
| *O. sativa* | `IRGSP-1.0_representative_2020-09-09.rev2.std_r_luc.gff3` | RAP-DB | **representative subset**, not the full annotation |
| *N. benthamiana* | `Niben101_annotation.gene_models.fix.std_r_luc.gff3` | SGN | scaffolds not chromosomes, one transcript per gene |
| *K. phaffii* | `GCF_000027005.1.ASM2700v1.chr.std_r_luc.gff3` | NCBI RefSeq | one transcript per gene |
| *S. cerevisiae* | **missing** | NCBI RefSeq | only `.gfd.gz` present, no GFF3 yet |

The *Chlamydomonas* file on the Mac is the v6.1 build rather than v6.0. The
two were compared and produce identical filtered counts, 16,883 genes and
31,858 transcripts, so results are comparable across machines. The chromosome
naming differs, `chr01` on WSL against `chromosome_01` on the Mac, which
matters for any `--seqids` filter.

### Still to obtain

- **GENCODE human**, needed for the validation gate and as an isoform-rich
  anchor. Annotation only, about 46 MB:
  `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz`
- **S. cerevisiae GFF3**, either exported from the existing `.gfd.gz` or
  downloaded from NCBI.
- **Full rice annotation** rather than the representative subset, if a source
  can be found.

## 7. Environment

| | Mac | WSL (`antec2025`) |
|---|---|---|
| conda environment | `ih` | `ibnu` |
| Python | 3.11.16 | 3.11 |
| hardware | Apple M5, 16 GB | Ryzen 9, 28 cores, 128 GB, 4 TB |
| role | development, writing, Level 0 | the parameter sweeps, later Level 1 |

Level 0 requires only Python 3.9 or newer and numpy, so it runs anywhere. The
heavier Level 1 stack, meaning an aligner and a quantifier, is intended for
the WSL machine, where Apptainer images already exist under
`/mnt/f/RA/Containers`.

On WSL, keep working data under `$HOME` rather than under `/mnt/`, which goes
through a translation layer and is noticeably slower for repeated reads.

## 8. Repository layout and the working rule

```
src/limen/   identify.py     the matrix, its rank, and diagnostic evidence
             annotation.py   GFF3 and GTF reading, standard library only
             toy.py          genes whose answers were derived by hand
test/        toy expectations and brute-force cross-checks
scripts/     batch jobs: produce data, write TSV, no plots
notebook/    interpret data, make figures, decide. reads what scripts wrote
ref/         manifest of genome sources, never the genomes themselves
docs/        METHODS.md, RESEARCH_PLAN.md
```

The division of labour is deliberate and worth preserving:

| | `scripts/*.py` | `notebook/*.ipynb` |
|---|---|---|
| does | produces data | interprets data |
| run | terminal, unattended, every genome | cell by cell, with you watching |
| input | annotation files | the TSVs the scripts wrote |
| output | TSV into `$LIMEN_WORK/results` | figures, tables, judgement |
| has | argparse, no plots, no hardcoded paths | plots, prints, exploration |
| heavy computation | yes | no, it loads results |

If it must be correct it belongs in `src/` with a test. If it must be looked
at it belongs in a notebook. If it is heavy enough to run unattended it
belongs in `scripts/` and writes a file the notebook then reads.

## 9. Running it

```bash
conda activate ibnu          # or ih on the Mac
cd ~/Code/limen
pip install -e ".[dev]"
pytest -q                    # expect 157 passing
```

One species:

```bash
python scripts/run_level0.py "$LIMEN_REF/<annotation>" \
    --species Chlamydomonas --jobs 12
```

All of them:

```bash
bash scripts/run_all.sh
```

The whole five-species panel takes about 30 seconds at twelve jobs. Each run
writes three files into `$LIMEN_WORK/results/level0`:

- `<species>.genes.tsv`, per gene and read length: rank, deficiency, whether
  identifiable, and the number of distinct read patterns
- `<species>.transcripts.tsv`, per transcript and read length: diagnostic
  start positions, total start positions, the fraction, and diagnostic reads
  per megabase of sequencing
- `<species>.run.json`, the sidecar recording the exact annotation file, its
  size, modification time and digest, every parameter, the parse statistics,
  timings and the package version

The sidecar exists so that any result can be traced back to the file and
parameters that produced it. A TSV without its sidecar should be treated as
unreliable.

## 10. Validation

Correctness is established in three ways, in increasing order of strength.
This matters more than usual here, because a simulation study that cannot
demonstrate its own correctness has no claim on anyone's attention.

**Toy genes with answers derived by hand.** The central case is four isoforms
produced by two independently skipped exons about 2000 nucleotides apart: rank
3 of 4 below 2002 nucleotides and rank 4 at or above it, where 2002 is the
shortest read that reaches from inside one skipped exon to inside the other.
A control gene has identical combinatorics with the exons 60 nucleotides
apart, and resolves at 62 nucleotides, which demonstrates that the effect is
about distance rather than isoform count.

**Brute-force cross-checks on randomised genes.** The fast segmentation
algorithms are compared against naive enumeration over thousands of random
gene structures, single-end and paired-end, on every test run.

**Published results.** Reproducing the human short-read identifiability rates
of Ferrer-Bonsoms and colleagues is treated as a prerequisite before the
read-length axis is extended beyond their range. It is a check on the
implementation, not a result. A second gate, reproducing the fragmentation
effect reported by Chen and colleagues, applies before any Level 1 claim.

### Bugs the tests have caught

Recorded because each produced plausible, wrong numbers that would not have
been noticed otherwise.

1. **Signatures are not constant between exon boundaries.** The first fast
   implementation assumed they were. Genomic coordinates shift with every
   base, so they are not. Thirty-seven of fifty-six tests failed immediately.
   The correct formulation turned out to be cleaner: what matters is not the
   signature but which transcripts share it, and that grouping genuinely is
   piecewise constant.
2. **Genomic segments could span an intron**, so read-start counts were
   overcounted. Fixed by making every exon edge a segment boundary.
3. **The transition point was bracketed on one side only.** A window starting
   exactly `R` before a boundary ends at the boundary without crossing it, so
   the change occurs at `b - R + 1` rather than `b - R`.
4. **The trailing run of start positions collapsed to width one**, because the
   final valid start had no boundary after it.
5. **Transcripts shorter than the read were counted by segment width** instead
   of once.
6. **Deduplication of identical transcript structures was nondeterministic.**
   It iterated a Python set, whose order varies between processes because of
   hash randomisation, so which transcript survived differed between runs.
   Caught by comparing serial against parallel output. For a project built on
   traceability this was the most damaging of the six.

### Things about the input that were not obvious

1. **Phytozome annotations contain no `exon` features at all.** Structure must
   be reconstructed from CDS and UTR features, which tile the mature
   transcript, merging pieces that abut. A parser reading only `exon` lines
   returns zero genes from these files, silently and without error.
2. **Annotation sources disagree about what a gene is.** NCBI files include
   `tRNA`, `ncRNA` and `pseudogene` entries, Araport includes transposons, and
   Phytozome marks `transposable_element_gene` separately. These are filtered
   deliberately rather than by default.
3. **An empty sweep is a result, not an error.** When no gene meets the
   isoform threshold the script writes empty tables and the sidecar and exits
   zero, so that a species with no annotated isoforms is recorded rather than
   aborting the run.

## 11. Prior work

The identifiability framework for short-read RNA sequencing is established:

- Hiller, D., Jiang, H., Xu, W., Wong, W.H. (2009) Identifiability of isoform
  deconvolution from junction arrays and RNA-Seq. *Bioinformatics* 25,
  3056-3059. Derives the conditions. Read length is not a variable.
- Ferrer-Bonsoms, J.A., Morales, X., Afshar, P.T., Wong, W.H., Rubio, A.
  (2022) On the identifiability of the isoform deconvolution problem.
  *Bioinformatics* 38, 1491-1496. Computes identifiability as a function of
  read and fragment length in order to select a library. Human only, 75 to 300
  nucleotides, unlimited depth assumed, long reads not considered. Code at
  `github.com/JFerrer-B/transcriptome-identifiability`, in R, covering
  chromosome 22 of GENCODE 24, with no licence stated.

Note that Wing Hung Wong is an author on both, so the line runs through one
group across thirteen years.

This repository contains an independent Python implementation written from the
published descriptions. No code is derived from those projects, which is a
requirement rather than a preference, since the reference repository states no
licence.

What is new here is the extension of the read-length axis into the long-read
range, the treatment of architecture as a variable, the replacement of binary
identifiability with a measure of how much evidence exists, and the
measurement of the distance between the structural ceiling and what real
pipelines recover.

The broader literature the project sits in is summarised with full citations
in `docs/RESEARCH_PLAN.md`. The papers that matter most are Chen et al. 2025
for the experimentally validated phantom isoforms, Han et al. 2024 for the
finding that short reads detect more splice junctions at matched coverage,
Apostolides et al. 2024 for joint short and long read quantification, and
Foord et al. 2023 for the statement of the underlying question about whether
features on a molecule are independent.

## 12. Open decisions

- Whether to analyse the current three usable species first, or obtain GENCODE
  and the full rice annotation before analysing. Analysing first would reveal
  within minutes whether gene-level architecture predicts anything at all,
  which is the assumption the whole reframe depends on.
- Whether to restrict every species to protein-coding genes on nuclear
  chromosomes, which would make them comparable but discard material.
- How to handle the organelle sequences present in Araport but absent from the
  *Chlamydomonas* chromosome-only file.
- Whether Level 1 generates reads with a purpose-built simulator in which
  every knob is an experimental variable, or wraps existing simulators whose
  platform presets bundle read length, error and artefacts together. The
  current intention is the former, with Badread used as a realism anchor
  rather than as a dependency.

## 13. Licence

MIT.
