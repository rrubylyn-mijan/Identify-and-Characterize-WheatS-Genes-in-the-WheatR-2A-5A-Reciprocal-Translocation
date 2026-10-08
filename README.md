# Identify-and-Characterize-WheatS-Genes-in-the-WheatR-2A-5A-Reciprocal-Translocation
The goal is to use Liftoff mappings to determine which WheatS  genes correspond to genomic regions involved in the WheatR reciprocal chromosome 2A/5A translocation, and then compare these genes with WheatG.

## At minimum, you need:
```bash
WheatS_genes_on_WheatR.gff3
WheatS_genes_on_WheatG.gff3
filtered_WheatS.gff3
```
## 1. Clean the names of gff3 file to match the genome naming
```bash
awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
$1 ~ /^chr[1-7][ABD]$/ {print}' \
TRAES.wheatR.pgsb.r1.Mar2024.high.gff3 \
> TRAES.wheatR.chromosomes.gff3

## Extract WheatS genes and wheatR genes
Make BED file:
cat > wheatS_chr5A-regions.bed <<EOF
chr5A	0	115000000	wheatS_chr5A_-0-115
EOF

cat > wheatS_chr2A-regions.bed <<EOF
chr2A	0	115000000	wheatS_chr2A-0-115
EOF

# do the same for wheatR

Extract genes:
ml bedtools2/2.31.1

bedtools intersect \
-a filtered_wheatS.gff3 \
-b wheatS_chr2A-regions.bed \
-wa \
> wheatS_chr2A_regions-0-115.gff3

## do the same for 5A and for wheatR

Get gene IDs:
grep -P "\tgene\t" wheatS_chr2A_regions-0-115.gff3 \
| cut -f9 \
| sed 's/.*ID=//' \
| sed 's/;.*//' \
> wheatS_chr2A_gene_ids-0-115.txt

## do the same for 5A and for wheatR
```

## 2. Get the coordinates of wheatS
```
grep -P "\tgene\t" wheatS_chr5A_regions-0-115.gff3 \
| awk 'BEGIN{OFS="\t"}{
split($9,a,";");
gsub("ID=","",a[1]);
print a[1],$1,$4,$5
}' \
> wheatS_chr5A-coordinates.tsv

## do the same for 5A. For wheatR, use the liftoff file but first generate it
```
## 3. Run liftoff to get the wheatS equivalent genes from wheatR translocated region
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

## 3. Use the saved gene list and the wheatR liftoff. 
```bash
# wheatS-on-wheatR-liftoff
grep -Ff wheatS_chr2A_gene_ids-0-115.txt wheatS_genes_on_wheatR.gff3 \
> wheatR-2A-0-115-liftoff-genes.gff3

# Extract the wheatS chr5A in liftoff
awk -F '\t' -v OFS='\t' '
BEGIN {
    print "gene_id", "chromosome", "start", "end"
}

# Read the gene IDs from the first file
FNR == NR {
    sub(/\r$/, "", $1)
    if ($1!="")wanted[$1] = 1
    next
}

# Read gene entries from the Liftoff GFF3
/^#/ || NF < 9 || $3 != "gene" { next }

{
    n = split($9, attributes, ";")
    id = ""

    for (i = 1; i <= n; i++) {
        if (attributes[i] ~ /^ID=/) {
            sub(/^ID=/, "", attributes[i])
            sub(/^gene:/, "", attributes[i])
            id = attributes[i]
            break
        }
    }

    if (id in wanted) {
        print id, $1, $4, $5
    }
}
' 5A-gene-IDs-sumai3 wheatS_genes_on_Rollag.gff3 \
> wheatS_5A_genes_mapped_to_wheatR_2A.tsv

## do the same for 5A
```

## 4. Make coordinate tables
```bash
# wheatR
grep -P "\tgene\t" wheatR-5A-0-115-liftoff-genes.gff3 \
| awk 'BEGIN{OFS="\t"}{
split($9,a,";");
gsub("ID=","",a[1]);
print a[1],$1,$4,$5
}' \
> wheatR-5A-0-115-coordinates.tsv

## do the same for 2A
```

## 5. Extract coding sequences (CDS) and translate them into proteins
```bash
ml miniconda3/25.5.1

conda init

source /apps/spack-managed-x86_64_v3-v1.1/gcc-11.5.0/miniconda3-25.5.1-75jzng7vis4hwdk3kzz5ywb4an56yzzz/etc/profile.d/conda.sh

conda activate gffread_env

# Clean the names of gff3 file to match the genome naming
awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
$1 ~ /^chr[1-7][ABD]$/ {print}' \
TRAES.wheatS.pgsb.r1.Mar2024.high.gff3 \
> TRAES.wheatS.chromosomes.gff3

## do the same for other wheat accessions

# Make wheatS proteins
gffread \
TRAES.wheatS.chromosomes.gff3 \
-g /directory/this/saved/wheatS_pm_v2.fasta \
-y wheats_proteins.fa

# Extract proteins:
ml seqkit/2.10.0

seqkit grep \
  -f <(awk '{print $1".1"}' wheats-gene-ids-5A-to-2A-0-115.txt) \
  wheats_proteins.fa \
  > wheats-gene-ids-5A-to-2A-0-115-proteins.fa
```

## 6. Run BLASTp
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p atlas
#SBATCH --mem=20GB
#SBATCH --time=72:00:00
#SBATCH -J blastp-gene-ids-5A-to-2A-0-115
#SBATCH -A genolabswheatphg

ml blastplus/2.17.0

cd /directory/this/saved/blastp-mapped-wheatr-sv

query="/directory/this/saved/blastp-mapped-wheatr-sv/wheats-gene-ids-5A-to-2A-0-115-proteins.fa"
out="wheats-gene-ids-5A-to-2A-0-115-blastp.tsv"

for attempt in 1 2 3 4 5
do
    echo "BLAST attempt $attempt"

    blastp \
    -query "$query" \
    -db swissprot \
    -remote \
    -out "$out" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
    -evalue 1e-10 \
    -max_target_seqs 5

    if [[ -s "$out" ]]; then
        echo "BLAST completed successfully"
        exit 0
    fi

    echo "BLAST failed or output is empty. Retrying in 5 minutes..."
    sleep 300
done

echo "BLAST failed after 5 attempts"
exit 1s

# Clean best-hit BLASTp result
sort -k1,1 -k6,6gr wheats-gene-ids-5A-to-2A-0-115-blastp.tsv | \
awk '!seen[$1]++' > wheats-gene-ids-5A-to-2A-0-115-function-besthit.tsv

# Make it easier to read
awk -F'\t' 'BEGIN{OFS="\t"; print "GeneID","SwissProt_hit","Percent_identity","Alignment_length","Evalue","Bitscore","BLAST_annotation"}
{
print $1,$2,$3,$4,$5,$6,$7
}' wheats-gene-ids-5A-to-2A-0-115-function-besthit.tsv > wheats-gene-ids-5A-to-2A-0-115-function-besthit-clean.tsv

# Remove transcript version so .1 and .2 become one gene ID
awk 'BEGIN{FS=OFS="\t"}
NR==1 {print; next}
{
  sub(/\.[0-9]+$/, "", $1)
  print
}' wheats-gene-ids-5A-to-2A-0-115-function-besthit-clean.tsv > wheats-gene-ids-5A-to-2A-0-115-besthit-clean-geneID.tsv

#  Keep only one BLAST hit per gene
sort -k1,1 -k6,6gr wheats-gene-ids-5A-to-2A-0-115-besthit-clean-geneID.tsv | \
awk 'BEGIN{FS=OFS="\t"} NR==1 || !seen[$1]++' \
> wheats-gene-ids-5A-to-2A-0-115-function-besthit-clean-geneID-unique.tsv

# Retain BLAST_annotation column for the genes listed in wheat-gene-ids-2A-to-5A-0-115.txt
awk -F'\t' 'BEGIN{OFS="\t"}
NR==FNR {
    if (NR>1) {
        match($0,/RecName: Full=([^;]+)/,m)
        if (m[1]!="")
            a[$1]=m[1]
    }
    next
}
{
    print $1, a[$1]
}' wheats-gene-ids-5A-to-2A-0-115-function-besthit-clean-geneID-unique.tsv wheats-gene-ids-5A-to-2A-0-115.txt \
> wheats-gene-ids-5A-to-2A-0-115-function-blastp-final-result.tsv


## do it for other wheat accessions you want to analyze
```

## 7. Run Interproscan
```bash
ml interproscan/5.78-109.0

interproscan.sh \
-i wheats-gene-ids-5A-to-2A-0-115-proteins.fa \
-f TSV \
-o wheats-gene-ids-5A-to-2A-0-115-interpro.tsv

# prepare interpro result
awk -F'\t' 'BEGIN{OFS="\t"}
NR==FNR {
    id=$1
    sub(/\.[0-9]+$/, "", id)

    # Only keep rows with an InterPro accession (IPRxxxxx)
    for(i=1;i<=NF;i++){
        if($i ~ /^IPR[0-9]+$/){
            ipr=$i
            desc=$(i+1)

            if(desc!="" && desc!="-"){
                if(interpro[id]=="")
                    interpro[id]=ipr": "desc
                else
                    interpro[id]=interpro[id]"; "ipr": "desc
            }
        }
    }
    next
}
{
    id=$1
    if(id in interpro)
        print id, interpro[id]
    else
        print id, "NA"
}' wheats-gene-ids-2A-to-5A-0-115-interpro.tsv wheats-gene-ids-2A-to-5A-0-115.txt \
> wheats-gene-ids-2A-to-5A-0-115-FINAL-InterPro.tsv

## do it with 5A to 2A
```

## 8. Get the log2foldchange
```bash
awk 'NR==FNR {if(NR>1)a[$1]=$2; next} NR==1 {print "Gene\tlog2FoldChange"; next} {print $1"\t"a[$1]}' logfold-wheats wheats-gene-ids-5A-to-2A-0-115.txt > wheats-5A-to-2A-0-115-s-g-r-expression-levels-trans-breakpoint-logfoldchange.tsv

## do the same for other wheat accessions
```

## 9. wheatR equivalent genes
```bash
# Save the gene IDs of wheatR
nano gene-IDs-wheat-2A-5A

wheatR chromosome       wheatR start    wheatR end
chr5A   251164  257767
chr5A   259218  263539
chr5A   279419  289395
chr5A   326283  327433
chr5A   389095  390999
chr5A   396052  400862
chr5A   471599  473099
chr5A   498871  501916
chr5A   511952  515630

# Run python script
python3 - <<'PY'
import csv
from collections import defaultdict

id_file = "wheats-gene-ids-5A-to-2A-0-115.txt"
liftoff_file = "wheatS_genes_on_wheatR.gff3"
rollag_file = "TRAES.WHEATR.chromosomes.gff3"
output_file = "wheatS_wheatR_ID_mapping-2A-5A.tsv"

def read_genes(path):
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 9 or fields[2] != "gene":
                continue

            attributes = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key] = value

            gene_id = attributes.get("ID")
            if gene_id:
                yield (
                    gene_id, fields[0],
                    int(fields[3]), int(fields[4]), fields[6]
                )

with open(id_file) as handle:
    requested = list(dict.fromkeys(
        line.strip() for line in handle if line.strip()
    ))

wanted = set(requested)
lifted = defaultdict(list)
native = defaultdict(list)

for gene in read_genes(liftoff_file):
    if gene[0] in wanted:
        lifted[gene[0]].append(gene)

for gene in read_genes(rollag_file):
    native[(gene[1], gene[4])].append(gene)

with open(output_file, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "Sumai3_gene_ID", "Rollag_gene_ID", "Status"
    ])

    for sumai_id in requested:
        if sumai_id not in lifted:
            writer.writerow([sumai_id, "NA", "NOT_FOUND_IN_LIFTOFF"])
            continue

        candidates = set()
        for _, chromosome, start, end, strand in lifted[sumai_id]:
            for rollag_id, _, rstart, rend, _ in native[(chromosome, strand)]:
                if max(start, rstart) <= min(end, rend):
                    candidates.add(rollag_id)

        if not candidates:
            writer.writerow([sumai_id, "NA", "NO_OVERLAPPING_ROLLAG_GENE"])
        else:
            status = (
                "SINGLE_OVERLAP_CANDIDATE"
                if len(candidates) == 1
                else "MULTIPLE_CANDIDATES_REVIEW"
            )
            for rollag_id in sorted(candidates):
                writer.writerow([sumai_id, rollag_id, status])

print(f"Saved: {output_file}")
print(f"Requested Sumai 3 genes: {len(requested)}")
print(f"Found in Liftoff: {len(lifted)}")
PY

```
