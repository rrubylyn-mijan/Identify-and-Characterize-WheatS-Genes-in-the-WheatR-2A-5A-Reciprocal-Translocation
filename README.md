# Identify-and-Characterize-WheatS-Genes-in-the-WheatR-2A-5A-Reciprocal-Translocation
The goal is to use Liftoff mappings to determine which WheatS  genes correspond to genomic regions involved in the WheatR reciprocal chromosome 2A/5A translocation, and then compare these genes with WheatG.

## At minimum, you need:
```bash
WheatS_genes_on_WheatR.gff3
WheatS_genes_on_Glenn.gff3
filtered_WheatS.gff3
```
## 1. Extract gene-level WheatS → WheatR mappings
```bash
awk -F'\t' '
BEGIN {OFS="\t"}

$3=="gene" {

    id=$9

    sub(/^.*ID=/,"",id)
    sub(/;.*/,"",id)

    print id,$1,$4,$5,$7
}' WheatS_genes_on_WheatR.gff3 \
> WheatS_genes_on_WheatR.coordinates.tsv

# Add a header:
sed -i \
'1iGeneID\tWheatR_chr\tWheatR_start\tWheatR_end\tStrand' \
WheatS_genes_on_WheatR.coordinates.tsv
```

## 2. From the result of your minimap alignment between WheatS and WheatR, identify the coordinates of translocated region. Look for large blocks where sequence from one chromosome aligns to a different chromosome
```bash
2A to 5A 251164 105638351
5A to 2A 30241 11232272
```

## 3. Extract the 2A → 5A translocated genes
```bash
awk -F'\t' -v OFS='\t' \
-v start=251164 \
-v end=105638351 '
NR==1 {next}

$2=="chr5A" &&
$3<=end &&
$4>=start &&
$1 ~ /\.2AG/ {

    print $1,$2,$3,$4,$5,"2A_to_5A"
}' WheatS_genes_on_WheatR.coordinates.tsv \
> WheatR_2A_to_5A_genes.tsv

# Add header:
sed -i \
'1iGeneID\tWheatR_chr\tWheatR_start\tWheatR_end\tStrand\tRegion' \
WheatR_2A_to_5A_genes.tsv
```

## 4. Extract the reciprocal 5A → 2A region
```bash
awk -F'\t' -v OFS='\t' \
-v start=30241 \
-v end=11232272 '
NR==1 {next}

$2=="chr2A" &&
$3<=end &&
$4>=start &&
$1 ~ /\.5AG/ {

    print $1,$2,$3,$4,$5,"5A_to_2A"
}' WheatS_genes_on_WheatR.coordinates.tsv \
> WheatR_5A_to_2A_genes.tsv

# Add header:
sed -i \
'1iGeneID\tWheatR_chr\tWheatR_start\tWheatR_end\tStrand\tRegion' \
WheatR_5A_to_2A_genes.tsv
```

## 5. Create a clean original WheatS GeneID
