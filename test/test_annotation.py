"""Loader behaviour, including reproducibility of transcript identity."""
import gzip

from limen.annotation import detect_format, load_annotation

GFF3 = """\
##gff-version 3
chr1\t.\tgene\t101\t2000\t.\t+\t.\tID=geneA
chr1\t.\tmRNA\t101\t2000\t.\t+\t.\tID=txA1;Parent=geneA
chr1\t.\texon\t101\t200\t.\t+\t.\tID=e1;Parent=txA1
chr1\t.\texon\t1801\t2000\t.\t+\t.\tID=e2;Parent=txA1
chr1\t.\tmRNA\t101\t2000\t.\t+\t.\tID=txA2;Parent=geneA
chr1\t.\texon\t101\t200\t.\t+\t.\tID=e3;Parent=txA2
chr1\t.\texon\t901\t1000\t.\t+\t.\tID=e4;Parent=txA2
chr1\t.\texon\t1801\t2000\t.\t+\t.\tID=e5;Parent=txA2
r_luc\t.\tgene\t1\t900\t.\t+\t.\tID=spikein
r_luc\t.\tmRNA\t1\t900\t.\t+\t.\tID=spiketx;Parent=spikein
r_luc\t.\texon\t1\t900\t.\t+\t.\tID=e6;Parent=spiketx
"""

PHYTOZOME = """\
##gff-version 3
chr1\t.\tgene\t101\t2000\t.\t+\t.\tID=G1
chr1\t.\tmRNA\t101\t2000\t.\t+\t.\tID=G1.1;Parent=G1
chr1\t.\tfive_prime_UTR\t101\t150\t.\t+\t.\tID=u1;Parent=G1.1
chr1\t.\tCDS\t151\t400\t.\t+\t0\tID=c1;Parent=G1.1
chr1\t.\tCDS\t1801\t2000\t.\t+\t0\tID=c2;Parent=G1.1
chr1\t.\tmRNA\t101\t2000\t.\t+\t.\tID=G1.2;Parent=G1
chr1\t.\tCDS\t101\t400\t.\t+\t0\tID=c3;Parent=G1.2
"""

# three transcripts, two of them structurally identical
DUPES = """\
##gff-version 3
chr1\t.\tgene\t1\t900\t.\t+\t.\tID=G
chr1\t.\tmRNA\t1\t900\t.\t+\t.\tID=G.zeta;Parent=G
chr1\t.\texon\t1\t100\t.\t+\t.\tID=z1;Parent=G.zeta
chr1\t.\texon\t801\t900\t.\t+\t.\tID=z2;Parent=G.zeta
chr1\t.\tmRNA\t1\t900\t.\t+\t.\tID=G.alpha;Parent=G
chr1\t.\texon\t1\t100\t.\t+\t.\tID=a1;Parent=G.alpha
chr1\t.\texon\t801\t900\t.\t+\t.\tID=a2;Parent=G.alpha
chr1\t.\tmRNA\t1\t900\t.\t+\t.\tID=G.beta;Parent=G
chr1\t.\texon\t1\t100\t.\t+\t.\tID=b1;Parent=G.beta
chr1\t.\texon\t401\t500\t.\t+\t.\tID=b2;Parent=G.beta
"""


def _write(tmp_path, name, text):
    p = tmp_path / name
    if name.endswith(".gz"):
        with gzip.open(p, "wt") as fh:
            fh.write(text)
    else:
        p.write_text(text)
    return str(p)


def test_gff3_with_exon_features(tmp_path):
    path = _write(tmp_path, "a.gff3", GFF3)
    assert detect_format(path) == "gff3"
    genes, st = load_annotation(path)
    assert len(genes) == 1 and genes[0].id == "geneA" and genes[0].n_iso == 2
    assert genes[0].transcripts[0].exons[0] == (100, 200)   # 1-based to 0-based
    assert st.skipped_excluded == 1                          # r_luc dropped
    assert st.tx_from_exon == 2 and st.tx_from_fallback == 0


def test_phytozome_cds_utr_fallback(tmp_path):
    """No exon features at all: structure must come from CDS and UTR, with
    abutting pieces merged into one exon."""
    path = _write(tmp_path, "p.gff3", PHYTOZOME)
    genes, st = load_annotation(path)
    assert st.tx_from_exon == 0 and st.tx_from_fallback == 2
    tx = {t.id: t for t in genes[0].transcripts}
    # 101-150 UTR abutting 151-400 CDS becomes a single exon
    assert tx["G1.1"].exons == ((100, 400), (1800, 2000))


def test_gzip_is_transparent(tmp_path):
    genes, _ = load_annotation(_write(tmp_path, "a.gff3.gz", GFF3))
    assert len(genes) == 1


def test_seqid_filter(tmp_path):
    genes, st = load_annotation(_write(tmp_path, "a.gff3", GFF3), seqids=["chr2"])
    assert genes == [] and st.skipped_seqid > 0


def test_exclude_can_be_disabled(tmp_path):
    genes, st = load_annotation(_write(tmp_path, "a.gff3", GFF3), exclude=[])
    assert st.skipped_excluded == 0
    assert any(g.seqid == "r_luc" for g in genes)


def test_duplicate_structures_collapse_deterministically(tmp_path):
    """Transcripts with identical exon chains are one observation unit, and
    which one represents them must not depend on hash ordering: a loader that
    returns different IDs on different runs makes results untraceable."""
    path = _write(tmp_path, "d.gff3", DUPES)
    genes, _ = load_annotation(path)
    assert len(genes) == 1
    ids = [t.id for t in genes[0].transcripts]
    assert len(ids) == 2, "the two identical structures should collapse to one"
    assert "G.alpha" in ids, "the lexicographically first ID should survive"
    assert "G.zeta" not in ids
    for _ in range(5):                       # stable across repeated loads
        assert [t.id for t in load_annotation(path)[0][0].transcripts] == ids
