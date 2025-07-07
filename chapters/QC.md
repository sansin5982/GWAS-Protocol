# Quality Control

**Quality control (QC)** in GWAS is a set of essential checks and
filters applied to both **samples and genetic variants** to remove
errors, technical artifacts, population biases, and unreliable data. QC
ensures that the statistical associations found are real, robust, and
replicable — not caused by poor genotyping, hidden relatedness,
population stratification, or data contamination. Without proper QC, a
GWAS can produce **false positives or miss true signals**, undermining
downstream analyses.

------------------------------------------------------------------------

## All required QC Steps

### Sample-Level QC:

This step is identify poorly or wrongly genotyped individuals in the
dataset. We want to identify bad samples before deciding which SNPs are
reliable. If we keep poor-quality samples in our dataset, they can
inflate missingness, cause false heterozygosity, or introduce spurious
HWE failures. Remove or fix **bad samples** → then we have a clean set
of individuals.

-   **Remove samples with high missing genotype rates**: Avoid
    poor-quality DNA or failed genotyping.

-   **Check and resolve sex mismatches**: Ensure genetic sex matches
    reported sex to detect sample swaps.

-   **Filter related or duplicate individuals**: Keep samples
    independent to prevent inflated test statistics.

-   **Identify ancestry outliers (PCA)**: Exclude individuals from
    different ancestral backgrounds to control population structure.

-   **Check autosomal heterozygosity**: Remove samples with unusually
    high or low heterozygosity, which may indicate contamination or
    inbreeding.

------------------------------------------------------------------------

### SNP-Level QC:

After removing bad samples, we try to remove SNPs wrongly genotyped or
have discrepency with an aim to find reliable SNPs. Once you know our
samples are good, we check each SNP for high missingness, low MAF, or
deviation from HWE. If we did SNP filters first, we’d be filtering based
on data that might still contain bad samples — so our SNP statistics
(missing rate, HWE) could be misleading.

-   **Filter SNPs with high missingness**: Remove poorly genotyped
    markers.

Apply Minor Allele Frequency (MAF) threshold: Exclude rare SNPs with
unstable estimates.

-   **Test Hardy-Weinberg Equilibrium (HWE) in controls**: Remove SNPs
    with deviations suggesting genotyping error.

-   **Filter poorly imputed SNPs**: Keep only well-imputed variants for
    reliable association tests.

------------------------------------------------------------------------

### Post-QC Check:

-   **Inspect QQ plots and genomic inflation factor (λGC)**: Confirm
    there is no hidden population stratification or systematic bias.
-   **Include principal components in association models**: Correct for
    residual population structure in the analysis.

------------------------------------------------------------------------

## Conclusion

QC is critical for valid, replicable GWAS results. Each step removes
noise, bias, or error, ensuring robust discovery and reliable downstream
biological interpretation.

------------------------------------------------------------------------

#### References

1- Marees, A.T., et al, 2018. A tutorial on conducting genome‐wide
association studies: Quality control and statistical analysis. *Int J
Methods Psychiatr Res*, Jun; 27(2): e1608.

2- Anderson, C.A. et al, 2010. Data quality control in genetic
case-control association studies. *Nat Protoc*, Sep:5(9):1564-73

3- Singh, Sandeep Kumar, “A Case-Only Genome-wide Association Study of
Gender- and Age-specific Risk Markers for Childhood Leukemia” (2015).
FIU Electronic Theses and Dissertations. 1832
