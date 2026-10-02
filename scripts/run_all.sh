#!/usr/bin/env bash
# Level 0 across the architecture panel. One TSV pair and sidecar per species.
set -euo pipefail
REF="${LIMEN_REF:?set LIMEN_REF}"
JOBS="${JOBS:-12}"

run () {  # species  filename
    echo "=============== $1"
    python scripts/run_level0.py "$REF/$2" --species "$1" --jobs "$JOBS"
}

run Komagataella  GCF_000027005.1.ASM2700v1.chr.std_r_luc.gff3
run Chlamydomonas CreinhardtiiCC_4532_707_v6.0.chr.std_r_luc.gff3
run Arabidopsis   Araport11_GFF3_genes_transposons.201606.std_r_luc.gff3
run Oryza         IRGSP-1.0_representative_2020-09-09.rev2.std_r_luc.gff3
run Nicotiana     Niben101_annotation.gene_models.fix.std_r_luc.gff3
