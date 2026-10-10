---
title: "Introduction to GWAS and Bioinformatics tools required"
output:
  md_document:
    variant: markdown
    preserve_yaml: true
---

# Introduction to GWAS and Bioinformatics tools required

A **genome-wide association study (GWAS)** is an observational genetic
study that tests whether genetic variants across the genome are
statistically associated with a defined phenotype. A **phenotype** is a
measurable characteristic, such as height, blood pressure, or whether a
person has a particular disease. A **genetic variant** is a difference
in DNA sequence among individuals. Many GWAS analyses focus on
**single-nucleotide polymorphisms (SNPs)**---positions at which
individuals may carry different DNA bases.

Rather than restricting attention to genes that researchers already
suspect, a GWAS systematically examines markers distributed throughout
the genome. For example, a study of blood pressure may detect an
association near a gene that was not previously considered part of
blood-pressure regulation. Such a result generates a biological
hypothesis for subsequent investigation; it does not demonstrate that
the nearest gene is responsible. This is often described as a
*genome-wide* or *hypothesis-generating* approach, but it is not
assumption-free: the phenotype, sampling strategy, statistical model,
ancestry adjustment, and quality-control decisions are all specified by
researchers.

**A central distinction:** GWAS identifies statistical *associations*;
it does not, by itself, establish *causation*. An associated marker may
influence the phenotype, may be correlated with a causal variant through
linkage disequilibrium (LD), or may reflect confounding or technical
artifacts. LD refers to the non-random association of alleles at
different genetic positions.

## From DNA variation to a phenotype

Our genome contains DNA sequences organized into chromosomes. At a given
SNP, an individual generally carries two alleles on an autosome, one
inherited from each parent. For example, if a SNP has alleles A and G,
the possible diploid genotypes are AA, AG, and GG. These genotypes may
be encoded as 0, 1, or 2 copies of a chosen effect allele in an additive
association model.

A genotype may contribute to a phenotype through biological pathways,
but phenotypes are also influenced by environment, age, sex, measurement
error, and many other genetic variants. The diagram below is a
conceptual illustration, **not evidence of a particular biological
pathway**.

::: {.figure style="text-align: center"}
`<img src="figures/ch01-fig-1-1-1.png" alt="Figure 1.1. Conceptual relationship between genetic variation and a phenotype. The arrows represent possible influences, not proof of causation."  />`{=html}
<p class="caption">
Figure 1.1. Conceptual relationship between genetic variation and a
phenotype. The arrows represent possible influences, not proof of
causation.
</p>
:::

## Why perform GWAS?

GWAS can identify genomic regions associated with disease susceptibility
or quantitative traits, prioritize hypotheses about biological
mechanisms, and generate summary statistics for subsequent research. For
instance, a variant associated with systolic blood pressure may point
researchers toward vascular, renal, or hormonal biology. Follow-up work
is needed to establish which mechanism, if any, explains the
association. Association results do **not** by themselves establish
clinical utility, predict disease accurately, or identify the causal
gene. Functional studies, replication, fine-mapping, and other evidence
may be needed.

## Candidate-gene studies versus GWAS

A **candidate-gene study** tests variants selected because of prior
biological expectations. A GWAS tests a much broader set of markers,
allowing discovery beyond known pathways. Both approaches can produce
false-positive results if study design or analysis is inadequate.

::: {.figure style="text-align: center"}
`<img src="figures/ch01-fig-1-2-1.png" alt="Figure 1.2. Candidate-gene studies test selected genomic regions; GWAS surveys markers across the genome. Marker density is schematic."  />`{=html}
<p class="caption">
Figure 1.2. Candidate-gene studies test selected genomic regions; GWAS
surveys markers across the genome. Marker density is schematic.
</p>
:::

# Study Design and Data Collection

## Define the phenotype before examining genotypes

A study needs a clear phenotype definition, inclusion and exclusion
criteria, an ascertainment strategy, and an analysis plan. Poor
phenotype classification can reduce power or bias an association
estimate. Our protocol distinguishes the following common designs:

  -----------------------------------------------------------------------
  Design            Outcome           Example           Common model
  ----------------- ----------------- ----------------- -----------------
  Case--control     Binary disease    Diabetes versus   Logistic
                    status            no diabetes       regression

  Quantitative      Continuous        Systolic blood    Linear regression
  trait             measurement       pressure          

  Time-to-event     Time until        Age at disease    Survival model,
                    event/censoring   onset             e.g. Cox
                                                        regression

  Case-only         Variation among   Early versus late Model depends on
  comparison        affected          onset among cases the case-only
                    individuals                         research question
  -----------------------------------------------------------------------

A case-only comparison does **not** estimate the same disease
susceptibility contrast as a population case--control study. Later
chapters will cover these designs in detail.

## A simple case--control example

Suppose we compare 100 cases with 100 controls at a SNP with alleles A
and G. For illustration, the counts below are fictional:

  Genotype     Cases   Controls
  ---------- ------- ----------
  AA              30         45
  AG              50         45
  GG              20         10
  Total          100        100

If G is the allele of interest, its frequency among cases is:

$$ f_{{G},\mathrm{cases}}=\frac{50+2(20)}{2(100)}=\frac{90}{200}=0.45.$$

Among controls:

$$ f_{{G},\mathrm{controls}}=\frac{45+2(10)}{2(100)}=\frac{65}{200}=0.325.$$

The G-allele frequency is therefore 45.0% among cases and 32.5% among
controls: a difference of **12.5 percentage points**. In everyday terms,
the G allele appears more frequently in the affected group in this small
fictional sample. That observation is a starting point for statistical
analysis, not a conclusion about disease risk. We must still ask whether
the difference is compatible with chance, whether population ancestry or
other factors could explain it, and how precisely an association model
estimates the effect. These counts are **illustrative**, not results
from our Hapmap data.

::: {.figure style="text-align: center"}
`<img src="figures/ch01-fig-1-3-1.png" alt="Figure 1.3. Genotype distribution in a fictional case–control study. The grouped bars show the number of cases and controls carrying AA, AG, or GG. These are illustrative counts, not Hapmap results."  />`{=html}
<p class="caption">
Figure 1.3. Genotype distribution in a fictional case--control study.
The grouped bars show the number of cases and controls carrying AA, AG,
or GG. These are illustrative counts, not Hapmap results.
</p>
:::

## The per-variant association test

For a binary phenotype, an additive logistic regression can be written
as:

$$\log\left(\frac{p}{1-p}\right)=\beta_{0}+\beta_{G} G+\beta_{A}\mathrm{Age}+\beta_{S}\mathrm{Sex}+\cdots$$

Here, $p$ is the modeled probability of being a case, $G$ is the number
of effect alleles (0, 1, or 2 for a diploid autosomal SNP), and
$\beta_G$ is the per-allele log-odds association coefficient. Its
exponent, $\exp(\beta_G)$, is the per-allele odds ratio under this
model. We will derive and interpret these quantities in the association
chapters; the formula is introduced here only to show what GWAS is
estimating.

For quantitative traits, an additive linear model commonly takes the
form $Y=\beta_0+\beta_GG+\cdots+\varepsilon$, where $\beta_G$ is the
estimated trait change per effect allele.

## Genome-wide testing and the Manhattan plot

A GWAS often tests hundreds of thousands or millions of markers. A
Manhattan plot displays each marker's genomic position on the horizontal
axis and $-\log_{10}(P)$ on the vertical axis. Smaller P values appear
higher. A commonly used threshold for common-variant GWAS is
$5\times10^{-8}$, although the appropriate threshold depends on the
analysis and testing burden. Statistical significance does not establish
causality.

The following figure uses **simulated P values solely for teaching**. It
is not an analysis of our HapMap dataset.

::: {.figure style="text-align: center"}
`<img src="figures/ch01-fig-1-4-1.png" alt="Figure 1.4. Illustrative Manhattan plot generated from simulated P values. The dashed line marks P = 5 × 10^-8; no biological discoveries are implied."  />`{=html}
<p class="caption">
Figure 1.4. Illustrative Manhattan plot generated from simulated P
values. The dashed line marks P = 5 × 10\^-8; no biological discoveries
are implied.
</p>
:::

## Data requirements

  -----------------------------------------------------------------------
  Input                   What it contains        Example
  ----------------------- ----------------------- -----------------------
  Genotypes               Alleles or genotype     PLINK `.bed/.bim/.fam`;
                          calls at variants       VCF

  Phenotypes              Disease status or       `.fam` phenotype column
                          measured trait          or separate phenotype
                                                  file

  Covariates              Variables used in       Age, sex, ancestry
                          adjustment              principal components,
                                                  batch where appropriate

  Metadata                Provenance and          Genome build, platform,
                          technical context       ancestry, sample
                                                  selection, file
                                                  versions
  -----------------------------------------------------------------------

**PLINK binary files:** `.bed` stores encoded genotypes, `.bim` stores
variant annotations including chromosome, identifier, position, and
allele codes, and `.fam` stores sample identifiers, pedigree/sex fields,
and a phenotype column. The two `.bim` allele columns should **not**
automatically be labeled REF and ALT. For a standard PLINK binary
case--control phenotype, `1` denotes control, `2` case, and `0` or `-9`
missing.

### Sample size and power

Statistical power depends on sample size, effect size, allele frequency,
phenotype prevalence and measurement quality, model choice, and
significance threshold. There is no universally sufficient sample size.
Small teaching datasets can demonstrate the workflow but are generally
underpowered for detecting the small effects common in complex-trait
GWAS.

### Population structure and confounding

Allele frequencies and phenotype distributions can differ across
ancestry groups for reasons unrelated to a causal effect at the tested
SNP. A poorly controlled comparison can therefore produce spurious
associations. We will study principal-component analysis (PCA), ancestry
outliers, relatedness, and covariate adjustment in their dedicated
chapters. Ancestry adjustment must be considered in relation to study
design; it is not a universal guarantee against confounding.

### Ethical considerations

Research using human genomic data requires appropriate consent, data
governance, privacy safeguards, and ethics review or a documented
exemption as applicable. Public teaching datasets do not eliminate the
need to respect licensing and data-use conditions.

# Material and Methods

## Our practical dataset: Hapmap

For simplicity, we will use the basename **`Hapmap`** throughout the
practical chapters. The files are derived from the European-ancestry
subset of HapMap Phase III genotypes distributed with the Marees et
al. (2018) tutorial. The original tutorial dataset is named
`HapMap_3_r3_1`; we have renamed the matching `.bed`, `.bim`, and `.fam`
files without changing their contents. Its binary phenotype is
**simulated**, not an observed disease outcome.

These are **processed genotype calls**, not raw SNP-array fluorescence
intensities. The dataset has already undergone processing before
distribution, so our QC is a reproducible *additional* QC workflow, not
a reconstruction of laboratory genotype calling.

Our recorded baseline PLINK log reports **165 people (80 males, 85
females), 1,457,897 variants, 112 nonmissing phenotype values, and an
overall genotyping rate of 0.997378**. The log also reports 225
heterozygous haploid calls, which require investigation in later
sex-chromosome QC. These figures describe the supplied tutorial files
and our initial inspection; no sample or variant has been excluded at
this stage.

The dataset contains autosomal and X-chromosome markers and
pseudoautosomal markers. It does not contain Y-chromosome markers coded
as chromosome 24 in the supplied `.bim`. Chromosome representation,
genome build, and allele conventions must be checked before any
downstream imputation.

## Directory layout

We keep the data and analysis outputs outside the GitHub Pages source
repository:

``` text
E:/Github_Classes/GWAS/
  01_Raw_Data/       # Starting Hapmap files (read-only baseline)
  02_Sample_QC/
  03_Variant_QC/
  04_Imputation/
  05_Association/
  06_Results/
  Logs/
  GWAS-Protocol/    # R Markdown chapters and publishable figures
```

The label `01_Raw_Data` means *raw input for our analysis*, not raw
laboratory measurements.

## Tools, software, and files required

**Required for the initial chapters:**

- **PLINK 1.9:** genotype inspection, QC, and introductory association
  analysis. Official site: <https://www.cog-genomics.org/plink/1.9/>.
- **R and RStudio (optional IDE):** summaries, figures, and reproducible
  R Markdown chapters. R: <https://www.r-project.org/>.
- **R Markdown/knitr:** rendering reproducible chapter text and plots.
- **Git:** tracking changes to the educational repository; large
  genotype files remain outside it.

**Introduced only when needed later:** PLINK 2, bcftools, samtools,
VCFtools, genome-build/liftover utilities, and imputation-preparation
tools. These are **not prerequisites** for Chapter 1 or for the initial
PLINK `.bed/.bim/.fam` inspection. Reference-panel and imputation-server
requirements can change and must be checked against the selected service
when we reach that chapter.

## Verifying the starting files in R

``` r
raw_dir <- "E:/Github_Classes/GWAS/01_Raw_Data"
files <- file.path(raw_dir, paste0("Hapmap", c(".bed", ".bim", ".fam")))
stopifnot(all(file.exists(files)))
file.info(files)[, c("size", "mtime")]

bim <- read.table(file.path(raw_dir, "Hapmap.bim"), header=FALSE)
fam <- read.table(file.path(raw_dir, "Hapmap.fam"), header=FALSE)
cat("Variants:", nrow(bim), "\n")
cat("Samples:", nrow(fam), "\n")
table(fam$V5, useNA="ifany")       # Recorded sex
table(fam$V6, useNA="ifany")      # Recorded phenotype
table(bim$V1, useNA="ifany")      # Chromosome codes
```

## Verifying the starting files with PLINK

From a terminal in the raw-data folder, the following command
**summarizes missingness without excluding variants or individuals**:

``` bash
./plink --bfile Hapmap --missing --out Raw_Missingness
```

It creates `Raw_Missingness.log`, `.imiss` (sample-level missingness),
and `.lmiss` (variant-level missingness). PLINK may also create `.hh` if
it detects heterozygous haploid calls. Our later QC chapters will record
the software version, exact commands, input/output paths, counts before
and after filters, warnings, decisions, and verification steps.

## The overall GWAS workflow

``` r
old <- par(mar=c(1,1,1,1));plot.new();plot.window(xlim=c(0,10),ylim=c(0,10))
labels <- c("Research question + phenotype definition",
            "Genotype acquisition + provenance",
            "Sample and variant QC",
            "Ancestry, relatedness + covariates",
            "Phasing/imputation (when applicable)",
            "Association testing + diagnostics",
            "Replication, interpretation + reporting")
ys <- seq(9.1,1.1,length.out=length(labels))
for (i in seq_along(labels)) {
 rect(1.4,ys[i]-.39,8.6,ys[i]+.39,col=if(i==2)"#FFF0CB" else "#D7E9F7",border="#315A78",lwd=1.25)
 text(5,ys[i],labels[i],cex=.96,font=if(i==2)2 else 1)
 if(i<length(labels)) arrows(5,ys[i]-.4,5,ys[i+1]+.4,length=.10,lwd=1.5)
}
```

::: {.figure style="text-align: center"}
`<img src="figures/ch01-fig-1-5-1.png" alt="Figure 1.5. Conceptual GWAS workflow. The Hapmap tutorial begins with already processed genotype calls, not raw laboratory data."  />`{=html}
<p class="caption">
Figure 1.5. Conceptual GWAS workflow. The Hapmap tutorial begins with
already processed genotype calls, not raw laboratory data.
</p>
:::

``` r
par(old)
```

The diagram presents the general order of work. Specific studies may
combine or adapt steps. Imputation is not compulsory for every
association analysis, and rigorous reporting continues throughout the
project.

## Reproducibility and troubleshooting record

For every computational step, we will record: (1) the input dataset and
genome build; (2) software and version; (3) the exact command or R code;
(4) input and output sample/variant counts; (5) thresholds and
rationale; (6) warnings and errors; (7) investigation and resolution;
and (8) post-step verification. This record will make it possible to
audit the practical analysis independently.

# Chapter summary

GWAS systematically tests genotype--phenotype associations across the
genome. Association does not by itself imply causation. The credibility
of a GWAS depends on phenotype quality, genotype QC, study design,
population structure, appropriate statistical models, multiple-testing
control, and transparent reporting. Our practical teaching dataset,
`Hapmap`, contains previously processed HapMap Phase III genotypes and a
simulated phenotype. We will use it to teach a reproducible PLINK-based
workflow, with all QC decisions documented.

# Concept checks

1.  Why can a statistically associated SNP fail to be causal?
2.  How does a case--control GWAS differ from a quantitative-trait GWAS?
3.  What is the difference between a SNP, an allele, and a genotype?
4.  Why is a `.bed/.bim/.fam` dataset not necessarily raw laboratory
    data?
5.  Why can population structure create a false association?
6.  What does a Manhattan plot show, and what does it *not* prove?
7.  Why do we retain the unmodified `Hapmap` starting files?

# References

1.  Marees AT, de Kluiver H, Stringer S, et al. A tutorial on conducting
    genome-wide association studies: Quality control and statistical
    analysis. *International Journal of Methods in Psychiatric
    Research.* 2018;27:e1608. <https://doi.org/10.1002/mpr.1608>.
2.  Visscher PM, Wray NR, Zhang Q, et al. 10 years of GWAS discovery:
    Biology, function, and translation. *American Journal of Human
    Genetics.* 2017;101:5--22.
    <https://doi.org/10.1016/j.ajhg.2017.06.005>.
3.  Abdellaoui A, Yengo L, Verweij KJH, Visscher PM. 15 years of GWAS
    discovery: Realizing the promise. *American Journal of Human
    Genetics.* 2023;110:179--194.
    <https://doi.org/10.1016/j.ajhg.2022.12.011>.
4.  National Human Genome Research Institute. Genome-Wide Association
    Studies Fact Sheet.
    <https://www.genome.gov/about-genomics/fact-sheets/Genome-Wide-Association-Studies-Fact-Sheet>.
5.  PLINK 1.9 documentation. <https://www.cog-genomics.org/plink/1.9/>.
