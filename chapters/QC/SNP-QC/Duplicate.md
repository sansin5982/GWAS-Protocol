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

### How can they break imputation?

This is critical:

Most phasing and imputation tools — like **SHAPEIT, EAGLE, Minimac4,
Beagle, IMPUTE2** — assume: \* **One unique site per genomic position**
\* Exactly **two alleles** at each position

**If a site has multiple alleles or repeats, the phasing step can’t
build a clean haplotype**.

#### Results:

-   Site might be silently skipped.
-   Haplotype model can break → segments fail to phase properly.
-   Reference panel matching fails → no imputed genotypes.
-   Or worse: the imputed data may contain strand flips or misaligned
    alleles.

Therefore: \* It’s strongly recommended to remove or fix these sites
before phasing and imputation — not after!

#### What should you do? Remove or keep?

Standard practice: \* For GWAS: **Remove duplicates and multi-allelics**
→ keep only unique bi-allelic SNPs. \* For phasing/imputation: **Must
remove or split multi-allelics** → tools expect bi-allelic input. \* For
post-imputation QC: Double-check the output for unresolved
multi-allelics if you merged multiple sources.

#### Check for duplicate sites

-   Duplicates break coordinate uniqueness — imputation tools fail or
    drop them.

<!-- -->

    ./plink --bfile Imiss_heter --list-duplicate-vars suppress-first --out Duplicate_SNPs

This will output a .dupvar file listing:

-   Duplicate IDs (same SNP name) or
-   Same position with multiple entries

#### Selecting only duplicated SNPs

    #### Bi allelic and Duplicate SNPs
    library(data.table)
    library(dplyr)

    # Load dupvar file
    dupvar <- fread("Duplicate_SNPs.dupvar") %>% 
      rename("SNP" = IDS, "BP" = POS)

    ./plink --bfile Imiss_heter --missing --out mydata_missing

This gives `mydata_missing.lmiss`

<table>
<thead>
<tr>
<th style="text-align: left;">CHR</th>
<th style="text-align: left;">SNP</th>
<th style="text-align: left;">N_MISS</th>
<th style="text-align: left;">N_GENO</th>
<th style="text-align: left;">F_MISS</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
</tr>
</tbody>
</table>

    # Load missingness data
    lmiss <- fread("mydata_missing.lmiss")
    bim <- fread("Imiss_heter.bim") %>% 
      rename("SNP" = V2, "BP" = V4) %>% 
      select(2, 4)

    dup_lmiss <- bim %>% 
      left_join(lmiss, by = "SNP") %>% 
      right_join(dupvar, by = "BP") %>% 
      select(-7, -9)


    # Group by BP and pick SNP with min F_MISS in each group
    dup_lmiss_best <- dup_lmiss[, .SD[which.min(F_MISS)], by = BP]


    to_exclude <- dup_lmiss_best %>% 
      select(CHR.x, BP, ALLELES, SNP.x) %>% 
      rename("CHR" = CHR.x, "POS" = BP, "IDS" = SNP.x)

    # Write to file
    write.table(
      to_exclude,
      "mydata_dropdups.txt",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )

**Best action**: Keep the best version (highest call rate) — or just
drop all duplicates if we don’t trust them.

#### Remove multi-allelic sites

Phasing/imputation assumes biallelic SNPs. Multi-allelic sites can’t be
modeled in the simple 0/1 allele framework.

     ./plink --bfile No_Dup_SNP --biallelic-only strict --make-bed --out biallelic

This step will create a new binary file containing only biallelic
alleles

#### References

Ricopili (PGC pipeline) — removes duplicates & non-biallelic sites
before per-SNP QC.

HRC Imputation SOP — same: bi-allelic only before imputation, QC must
match.

1000 Genomes and UK Biobank pipelines: clean VCFs must be bi-allelic for
phasing & all QC steps.
