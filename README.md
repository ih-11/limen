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
