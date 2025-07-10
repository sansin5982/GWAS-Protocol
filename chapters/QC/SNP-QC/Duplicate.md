# Duplicate, Multi-allelic, and Non-biallelic sites

### Duplicate sites

-   Two or more SNPs that share the **same chromosome and base
    position**.

-   Example: Our dataset has `rs123` and `rs456` both mapped to
    `Chr 1:1234567`.

-   These might come from:

    -   Genotyping array annotation overlaps
    -   Different strand versions of same site
    -   Bad merges or duplicate records

### Multi-allelic or non-biallelic sites

-   **Bi-allelic** means a SNP has **one reference allele** and **one
    alternate allele** (e.g., A/G).
-   **Multi-allelic** means there are **more than two alleles at the
    same genomic position**.
    -   Example: at Chr 1:1234567 the possible alleles are A, G, and T.
-   Indels (insertions/deletions) or complex variants can also introduce
    multi-allelic calls.
-   In VCFs, these show up with multiple ALT fields, like ALT=G,T.

#### Key point:

Modern **GWAS tools, QC tools, phasing tools, and imputation servers
assume each site is unique and bi-allelic** — by default.

## Why do they create problems for GWAS?

#### Technical issues:

-   Many GWAS tools (e.g., PLINK) assume SNPs are **bi-allelic** when
    coding genotypes.
-   If we have more than two alleles, PLINK must decide how to encode
    them as 0, 1, 2 — but this is ambiguous for &gt;2 alleles.
-   Some older GWAS software can crash, error out, or silently miscode
    genotypes.

#### Statistical issues:

-   Association tests assume our variant is well-defined:

    -   The genotypes map cleanly to **homozygous reference,
        heterozygous**, and **homozygous alternate**.

-   Multi-allelic SNPs can split our sample into multiple genotype
    groups with very small counts → unstable p-values.

-   Duplicate positions may inflate the apparent number of tests →
    messes up our multiple testing correction.

#### Downstream QC confusion:

-   Hardy-Weinberg test, allele frequency checks, and relatedness stats
    can get distorted by multi-allelics and duplicated markers.

#### References

Ricopili (PGC pipeline) — removes duplicates & non-biallelic sites
before per-SNP QC.

HRC Imputation SOP — same: bi-allelic only before imputation, QC must
match.

1000 Genomes and UK Biobank pipelines: clean VCFs must be bi-allelic for
phasing & all QC steps.
