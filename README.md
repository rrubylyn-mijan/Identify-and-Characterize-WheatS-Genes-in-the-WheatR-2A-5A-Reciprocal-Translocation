# Identify-and-Characterize-WheatS-Genes-in-the-WheatR-2A-5A-Reciprocal-Translocation
The goal is to use Liftoff mappings to determine which WheatS  genes correspond to genomic regions involved in the WheatR reciprocal chromosome 2A/5A translocation, and then compare these genes with WheatG.

## At minimum, you need:
```bash
WheatS_genes_on_WheatR.gff3
WheatS_genes_on_Glenn.gff3
filtered_WheatS.gff3
```
## 1. Extract WheatS genes from translocated region
```bash
Make BED file:
cat > wheatR_translocation_regions.bed <<EOF
chr5A	0	115000000	wheatR_chr5A_to_wheatR_chr2A-0-115
chr2A	0	115000000	wheatR_chr2A_to_wheatR_chr5A-0-115
EOF

Extract genes:
bedtools intersect \
-a filtered_wheatS.gff3 \
-b wheatR_translocation_regions.bed \
-wa \
> wheatS_wheatR_translocation_regions-0-115.gff3

Get gene IDs:
grep -P "\tgene\t" wheatS_wheatR_translocation_regions-0-115.gff3 \
| cut -f9 \
| sed 's/.*ID=//' \
| sed 's/;.*//' \
> wheatS_wheatR_translocation_gene_ids-0-115.txt
```

## 2. Run liftoff to get the wheatS equivalent genes from wheatR translocated region
```bash
nano wheatS_genes_on_wheatR.sh
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH -p atlas
#SBATCH --mem=100GB
#SBATCH --time=24:00:00
#SBATCH -A genolabswheatphg

# ============================================================
# LIFTOFF: wheatS -> wheatR
#
# Purpose:
# Transfer the wheatR high-confidence gene annotation
# onto the wheatR genome assembly.
#
# SOURCE:
#   wheatS genome
#   wheatS PGSB high-confidence GFF3
#
# TARGET:
#   wheatR genome
#
# Main output:
#   wheatS_genes_on_wheatR.gff3
# ============================================================


# ============================================================
# 1. STOP IF ANY COMMAND FAILS
# ============================================================

set -euo pipefail


# ============================================================
# 2. LOAD SOFTWARE
# ============================================================

ml miniconda3/25.5.1

source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate liftoff

ml minimap2/2.30


# ============================================================
# 3. CHECK SOFTWARE
# ============================================================

echo "============================================================"
echo "SOFTWARE"
echo "============================================================"

echo "Liftoff:"
which liftoff
liftoff --version

echo

echo "Minimap2:"
which minimap2
minimap2 --version

echo "============================================================"
echo


# ============================================================
# 4. INPUT FILES
# ============================================================

# wheatS genome assembly
wheatS_FA="/directory/this/saved/wheatS_pm_v2.fasta"

# wheatS PGSB high-confidence annotation
wheatS_GFF="/directory/this/saved/3_gff_extracts/filtered_wheatS.gff3"

# wheatR genome
wheatR_FA="/directory/this/saved/wheatR_pm_v1.fasta"


# ============================================================
# 5. OUTPUT DIRECTORY
# ============================================================

OUT_DIR="/directory/this/saved/1lift-off-wheatS-wheatR"

mkdir -p "$OUT_DIR"


# ============================================================
# 6. LOCAL SCRATCH DIRECTORY
# ============================================================

LOCAL_WORK="$TMPDIR/liftoff_wheatS_to_wheatR"

mkdir -p "$LOCAL_WORK"


# ============================================================
# 7. PRINT JOB INFORMATION
# ============================================================

echo "============================================================"
echo "LIFTOFF: wheatS -> wheatR"
echo "============================================================"

echo "SLURM Job ID:       ${SLURM_JOB_ID:-NA}"
echo "CPUs:               ${SLURM_CPUS_PER_TASK:-1}"
echo

echo "SOURCE FASTA:"
echo "$wheatS_FA"
echo

echo "SOURCE GFF3:"
echo "$wheatS_GFF"
echo

echo "TARGET FASTA:"
echo "$wheatR_FA"
echo

echo "LOCAL WORK:"
echo "$LOCAL_WORK"
echo

echo "FINAL OUTPUT:"
echo "$OUT_DIR"

echo "============================================================"
echo


# ============================================================
# 8. CHECK INPUT FILES
# ============================================================

echo "Checking input files..."

for FILE in "$wheatS_FA" "$wheatS_GFF" "$wheatR_FA"; do

    if [[ ! -s "$FILE" ]]; then
        echo
        echo "ERROR: Input file is missing or empty:"
        echo "$FILE"
        exit 1
    fi

    echo "FOUND: $FILE"

done

echo
echo "All input files found."
echo


# ============================================================
# 9. COUNT SOURCE GENES
# ============================================================

N_SOURCE_GENES=$(awk -F'\t' \
    '$3=="gene"{n++} END{print n+0}' \
    "$wheatS_GFF")

echo "Number of genes in wheatS source GFF3:"
echo "$N_SOURCE_GENES"

echo


# ============================================================
# 10. RUN LIFTOFF
# ============================================================

echo "============================================================"
echo "STARTING LIFTOFF"
echo "============================================================"

date

echo

liftoff \
    -g "$wheatS_GFF" \
    -o "$LOCAL_WORK/wheatS_genes_on_wheatR.gff3" \
    -u "$LOCAL_WORK/wheatS_unmapped_features.txt" \
    -dir "$LOCAL_WORK/intermediate_files" \
    -p "$SLURM_CPUS_PER_TASK" \
    "$wheatR_FA" \
    "$wheatS_FA"

echo
date

echo "============================================================"
echo "LIFTOFF COMMAND FINISHED"
echo "============================================================"
echo


# ============================================================
# 11. VERIFY MAIN OUTPUT
# ============================================================

if [[ ! -s "$LOCAL_WORK/wheatS_genes_on_wheatR.gff3" ]]; then

    echo "ERROR:"
    echo "Liftoff did not produce a non-empty GFF3:"
    echo "$LOCAL_WORK/wheatS_genes_on_wheatR.gff3"

    exit 1

fi

echo "Mapped GFF3 successfully produced."
echo


# ============================================================
# 12. COPY MAIN OUTPUT TO PERMANENT STORAGE
# ============================================================

cp \
    "$LOCAL_WORK/wheatS_genes_on_wheatR.gff3" \
    "$OUT_DIR/wheatS_genes_on_wheatR.gff3"

echo "Copied mapped GFF3 to:"
echo "$OUT_DIR/wheatS_genes_on_wheatR.gff3"

echo


# ============================================================
# 13. COPY UNMAPPED FILE IF IT EXISTS
# ============================================================

if [[ -s "$LOCAL_WORK/wheatS_unmapped_features.txt" ]]; then

    cp \
        "$LOCAL_WORK/wheatS_unmapped_features.txt" \
        "$OUT_DIR/wheatS_unmapped_features.txt"

    echo "Copied unmapped features file to:"
    echo "$OUT_DIR/Sumai3_unmapped_features.txt"

else

    echo "No non-empty unmapped-features file was produced."

fi

echo


# ============================================================
# 14. BASIC OUTPUT QC
# ============================================================

MAPPED_GFF="$OUT_DIR/wheatS_genes_on_wheatR.gff3"

N_MAPPED_GENES=$(awk -F'\t' \
    '$3=="gene"{n++} END{print n+0}' \
    "$MAPPED_GFF")


echo "============================================================"
echo "LIFTOFF SUMMARY"
echo "============================================================"

echo "Source wheatS genes:     $N_SOURCE_GENES"
echo "Mapped gene records:      $N_MAPPED_GENES"

echo
echo "Mapped GFF3:"
echo "$MAPPED_GFF"

if [[ -s "$OUT_DIR/wheatS_unmapped_features.txt" ]]; then

    N_UNMAPPED=$(wc -l < \
        "$OUT_DIR/wheatS_unmapped_features.txt")

    echo
    echo "Unmapped feature entries: $N_UNMAPPED"
    echo
    echo "Unmapped features:"
    echo "$OUT_DIR/wheatS_unmapped_features.txt"

fi

echo
echo "============================================================"
echo "JOB COMPLETED SUCCESSFULLY"
echo "============================================================"

date
```

## 3. Use the saved gene list and the wheatR liftoff
```bash
# wheatS-on-wheatR-liftoff
grep -Ff wheatR-gene-ids-5A-to-2A-0-115.txt wheatS_genes_on_wheatR.gff3 \
> wheatR-5A-to-2A-0-115-liftoff-genes.gff3

# wheatS-equivalent
grep -Ff wheatR-gene-ids-5A-to-2A-0-115.txt filtered_wheatS.gff3 \
> wheatS-equivalent-5A-to-2A-0-115-liftoff-genes.gff3
```

## 4. Make coordinate tables
```bash
# wheatR
grep -P "\tgene\t" wheatR-5A-to-2A-0-115-liftoff-genes.gff3 \
| awk 'BEGIN{OFS="\t"}{
split($9,a,";");
gsub("ID=","",a[1]);
print a[1],$1,$4,$5
}' \
> wheatR-5A-to-2A-0-115-coordinates.tsv

# wheatS-equivalent
grep -P "\tgene\t" wheatS-equivalent-5A-to-2A-0-115-liftoff-genes.gff3 \
| awk 'BEGIN{OFS="\t"}{
split($9,a,";");
gsub("ID=","",a[1]);
print a[1],$1,$4,$5
}' \
> wheatS-equivalent-5A-to-2A-0-115-coordinates.tsv
```

## 5. Create a clean original WheatS GeneID
