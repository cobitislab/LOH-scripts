# LOH-scripts

Collection of R & python scripts used in our recent project exploring the Loss of Heterozygosity (LOH) in clonal *Cobitis* fish hybrids.

**Preprint** on [bioRxiv](https://doi.org/10.1101/2020.07.30.229369).

## Pairwise mismatches
The script `pairwise_mismatches.r` counts mismatches in genotypes between pairs of samples.

**Input:** The script is designed to work with output of GATK `VariantsToTable` tool. Minimal input table requires three types of fields - `CHROM`, `POS`, and `GT`. The table can have other columns (such as `DP`, `AD`, etc.) - the script only works with `GT` columns and ignores the rest.

Example command to obtain such table:
```bash
gatk3 -T VariantsToTable -R ref.fa -V var.vcf -F CHROM -F POS -GF GT -o table.tsv
```
