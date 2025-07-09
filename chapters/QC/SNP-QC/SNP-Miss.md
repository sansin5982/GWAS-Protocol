# SNP Call rate

SNP call rate means:

-   For each SNP, what **percentage of samples** have a valid genotype
    call.
-   The flip side is **missingness** — how many samples failed to get a
    genotype call for that SNP.

#### Example:

We have 1,000 samples.

-   SNP\_A is successfully genotyped in 990 samples → Call rate = 99% →
    Missingness = 1%.

-   SNP\_B has only 940 good calls → Call rate = 94% → Missingness = 6%.

#### Why filter by SNP call rate?

SNPs with **lots of missing calls** often fail due to:

-   Poor probe performance on the genotyping chip.
-   Lab contamination.
-   Batch-specific technical failures.

Including them can:

-   Add noise.

-   Bias allele frequency estimates.

-   Cause false positives if missingness correlates with case/control
    status.

**Good genotyping labs usually produce SNPs with &gt;98–99% call
rates.**

-   So typical GWAS filters remove SNPs missing in more than 2–5% of
    samples.

# References

1- Marees, A.T., et al, 2018. A tutorial on conducting genome‐wide
association studies: Quality control and statistical analysis. *Int J
Methods Psychiatr Res*, Jun; 27(2): e1608.

2- Anderson, C.A. et al, 2010. Data quality control in genetic
case-control association studies. *Nat Protoc*, Sep:5(9):1564-73

3- Singh, Sandeep Kumar, “A Case-Only Genome-wide Association Study of
Gender- and Age-specific Risk Markers for Childhood Leukemia” (2015).
FIU Electronic Theses and Dissertations. 1832
