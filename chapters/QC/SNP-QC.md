# Per SNP quality control

**Per-SNP quality control** means evaluating the **quality and
reliability** of each genetic variant (SNP) in our dataset before
testing it for association with our trait of interest.

In any GWAS, we may start with **hundreds of thousands to millions** of
SNPs across the genome. But not every SNP is equally reliable — some may
have been poorly genotyped, rarely observed, or show patterns that
suggest technical errors.

So, per-SNP QC systematically **filters out bad or suspicious
variants**, leaving only high-confidence markers for robust downstream
analyses.

### Why is Per-SNP QC Necessary?

#### 1. Avoid Technical Errors:

Some SNPs consistently have poor call rates or high missingness.
Including these can introduce **genotyping noise**, increasing false
positives or masking true signals.

#### 2. Exclude Rare, Unreliable Signals:

Extremely rare variants in small samples often lack enough data for
stable statistical testing. They can inflate standard errors or lead to
unreliable results if left unfiltered.

#### 3. Detect Systematic Problems:

A SNP whose genotype frequencies **deviate too much** from what’s
expected (like Hardy-Weinberg Equilibrium in controls) may indicate:

-   Lab contamination.
-   Genotyping assay problems.
-   DNA sample mix-ups.
-   Removing these ensures only **biologically plausible** variants
    remain.

#### 4. Protect Downstream Analyses:

GWAS assumes that each tested SNP is a valid, well-measured marker.
Poor-quality SNPs distort:

-   Allele frequencies,
-   Effect estimates,
-   P-values, and can lead to irreproducible results.

Per-SNP QC is about trust:

> Are we confident this SNP really shows what we think it does?

A robust GWAS filters out unreliable SNPs **before** the main analysis —
ensuring that any significant hits are due to **real biological
effects**, not technical artifacts.

Per-SNP QC is a fundamental safeguard in modern genetic studies. Without
it, even the most advanced association methods can produce misleading
results. It’s simple in principle: keep only high-quality, well-behaved
SNPs, so our study has the best chance of finding true genetic
associations that replicate.

#### SNP call rate (Missingness)

Remove SNPs with too much missing data across all samples. High
missingness suggests unreliable genotyping at that position. Remove SNPs
missing in &gt;2–5% of samples.

#### Minor Allele Frequency

Filter out SNPs that are too rare to analyze robustly. Variants with
extremely low frequency often have unstable allele frequency estimates
and low statistical power. Drop SNPs with MAF below 0.01 (1%) or 0.05
(5%) — this depends on our sample size and study design.

#### SNPs failing Hardy Weinberg equilibrium

Test whether observed genotype frequencies match those expected under
random mating. Large deviations in controls often signal genotyping
errors, population substructure, or sample contamination. Remove SNPs
with HWE p-value below a strict threshold (e.g., 1e-6) in controls only
— not in cases, because true disease associations can also cause
deviations.

#### Duplicates, Multi-allelic or Non-biallelic Sites

Many QC workflows keep only biallelic SNPs. Most standard GWAS tools
assume SNPs have only two alleles (e.g., A/T). Complex variants often
require special treatment or separate analysis.

Though we can remove all low-quality SNPs step by step, but filtering
them together in a single PLINK run is standard, efficient, and
recommended.

#### The filters are independent

-   Each QC filter (SNP call rate, MAF, HWE) is just a logical check.
-   If a SNP fails any one filter, you want it gone.

#### One step is faster and easier to track

-   Running everything together means:

    -   Only one PLINK command.
    -   One consistent output.
    -   Less chance of accidental reprocessing or forgetting an
        intermediate file.

#### Exception: When you must filter on a subset

-   One case when we DO break it up:

-   **HWE**: We should filter on controls only, not all samples.

So we might:

-   Extract controls → run `--hwe` → get list of SNPs that pass.⃣
-   Remove those SNPs from the full dataset.

# References

1- Marees, A.T., et al, 2018. A tutorial on conducting genome‐wide
association studies: Quality control and statistical analysis. *Int J
Methods Psychiatr Res*, Jun; 27(2): e1608.

2- Anderson, C.A. et al, 2010. Data quality control in genetic
case-control association studies. *Nat Protoc*, Sep:5(9):1564-73

3- Singh, Sandeep Kumar, “A Case-Only Genome-wide Association Study of
Gender- and Age-specific Risk Markers for Childhood Leukemia” (2015).
FIU Electronic Theses and Dissertations. 1832

4- Ricopili pipeline (PGC): does HWE control-only first → then combines
other filters.
