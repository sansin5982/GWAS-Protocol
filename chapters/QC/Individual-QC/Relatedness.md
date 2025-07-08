<script type="text/javascript" async
    src="https://polyfill.io/v3/polyfill.min.js?features=es6">
</script>
<script type="text/javascript" async
    src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.0/es5/tex-mml-chtml.js">
</script>

# Relatedness

**Relatedness QC** is the step in genome-wide association studies where
we **identify and manage individuals who are closely related to each
other** in our sample.

Most GWAS assume that all individuals are **statistically independent**
— meaning they’re **unrelated** beyond distant population-level
relationships. If we leave related samples in the dataset, they can
cause **inflated test statistics**, spurious associations, or biased
heritability estimates.

### Why Related Individuals Cause Problems

-   **Non-independence**: Related individuals share long stretches of
    DNA identical by descent (IBD). This means their genotypes are not
    fully independent observations.

-   **Allele frequency bias**: The presence of close relatives can
    over-represent certain alleles in our sample.

-   **Inflated test statistics**: The GWAS model assumes that the errors
    are independent. Relatedness violates this, so our p-values may be
    **too small** — leading to false positives.

#### How is Relatedness Measured?

The **relatedness coefficient (or pi-hat)** estimates the proportion of
the genome that two individuals share **identical by descent (IBD)**.

**Standard thresholds**:

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr>
<th>PIHAT values</th>
<th>Relation</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>Indicates identical twins or duplicates, where individuals share
100% of their genome IBD.</td>
</tr>
<tr>
<td>0.5</td>
<td>Suggests first-degree relatives, like parent and child or siblings,
sharing approximately 50% of their genome IBD.</td>
</tr>
<tr>
<td>0.25</td>
<td>Indicates second-degree relatives, like grandparents and
grandchildren or half-siblings, sharing about 25% of their genome
IBD.</td>
</tr>
<tr>
<td>0.125</td>
<td>Indicates third-degree relatives, like first cousins, sharing
roughly 12.5% of their genome IBD</td>
</tr>
</tbody>
</table>

### Identical By Descent (IBD)

#### Definition:

Two alleles are **Identical By Descent (IBD)** if they are **exact
copies inherited from the same ancestor without any mutation in
between**.

#### Key point:

IBD means there is a **known common ancestor** in the family tree — the
alleles are **genetically “linked” through inheritance**.

#### Example:

If you and your sibling both have an “A” allele at a locus, and you both
got it from your mother, that “A” is IBD.

### When is IBD used?

-   To measure **relatedness** in families: siblings, parent–child,
    cousins.
-   To detect **pedigree structure**.
-   In **pi-hat calculations**, where you want to estimate how much DNA
    two people share because of shared ancestry.
-   For **linkage studies**, where you track inheritance in families to
    find disease genes.

<img src="IBD.png" alt="Identical By Descent" width="1500" />
<p class="caption">
Identical By Descent
</p>

[Image
Credit](https://bio.libretexts.org/Bookshelves/Genetics/Population_and_Quantitative_Genetics_%28Coop%29/02%3A_Allele_and_Genotype_Frequencies)

### Identical By State (IBS)

Two alleles are **Identical By State (IBS)** if they are **the same
allele type**, regardless of whether they come from the **same
ancestor**.

IBS only looks at the current state of the allele — not how it got
there.

#### Example:

A random person ad I both have an “A” at a SNP — but maybe I got mne
from my mom’s lineage and they got theirs from an unrelated lineage.
**These two alleles are IBS, but not necessarily IBD**.

#### When is IBS used?

-   To compute **genetic similarity** between any two individuals, even
    if they’re unrelated.
-   IBS is the **observed data** — the raw count of matching alleles.
-   Tools like PLINK compare IBS at each SNP, then use population allele
    frequencies to infer IBD.

#### How They Relate

-   **IBS is what we observe directly**: “Do these two individuals have
    the same allele?”
-   **IBD is what we infer statistically**: “Are they the same because
    they both inherited it from a common ancestor?”

> So, all IBD alleles are also IBS — but not all IBS alleles are IBD.

#### Analogy:

Imagine two people wearing the same shirt.

-   If they got it from the **same parent** who handed it down — that’s
    **IBD**.
-   If they both just happen to buy the same shirt in a store — that’s
    **IBS**.

#### How They’re Estimated

-   Tools like PLINK, KING, and GCTA compare genotypes across thousands
    of SNPs.
-   They start by calculating IBS — simple allele matches.
-   Then they **adjust for population allele frequencies** to estimate
    whether the matching is just random (IBS) or likely due to
    inheritance (IBD).

#### Why the Distinction Matters

-   In **GWAS QC**, we care about **IBD** — we want to know who’s
    actually related.
-   IBS alone can’t distinguish between actual relatives and unrelated
    people with coincidentally similar genotypes.
-   Correct IBD estimates ensure you don’t accidentally include
    relatives in an “unrelated” study, which would **inflate association
    results**.

#### Summary

<table>
<colgroup>
<col style="width: 13%" />
<col style="width: 45%" />
<col style="width: 40%" />
</colgroup>
<thead>
<tr>
<th>Feature</th>
<th>IBD</th>
<th>IBS</th>
</tr>
</thead>
<tbody>
<tr>
<td>Means?</td>
<td>Same allele from same ancestor</td>
<td>Same allele, origin irrelevant</td>
</tr>
<tr>
<td>Connection</td>
<td>Inherited through a known pedigree</td>
<td>Could be chance or related</td>
</tr>
<tr>
<td>Used for</td>
<td>Pedigrees, relatedness, linkage</td>
<td>Raw match count, similarity</td>
</tr>
<tr>
<td>Example</td>
<td>Shared segment from grandparent</td>
<td>Same SNP allele in population</td>
</tr>
</tbody>
</table>

#### Applications of Pi-Hat Values

The **pi-hat** value, or **pairwise IBD (Identity by Descent)
estimate**, is a widely used measure that quantifies how much of the
genome two individuals share because they inherited it from a common
ancestor. It has multiple practical applications in genetic analysis,
quality control, and population research.

#### 1. Pedigree Reconstruction

Pi-hat values help researchers **reconstruct family trees** and infer
unknown or partial relationships among individuals in genetic data.

#### How it works:

By comparing pi-hat estimates, you can categorize pairs into known
familial classes:

-   ~0.5 → parent-offspring or full siblings
-   ~0.25 → half-siblings, grandparent-grandchild,
    aunt/uncle-niece/nephew
-   ~0.125 → first cousins

If we have no recorded pedigree or incomplete family records, pi-hat
helps identify relatives in a population-based sample. This is
especially useful in large biobank datasets, ancient DNA studies, or
conservation genetics where family links are missing or uncertain.

#### 2. Genome-Wide Association Studies (GWAS)

In GWAS, pi-hat is essential for **detecting related individuals** whose
genetic similarity can violate the assumption of **sample
independence**.

#### How it works:

If close relatives remain in the data, they can **inflate test
statistics**, causing false positive associations or misleading
p-values. By identifying pairs with high pi-hat values, researchers can:

-   Remove one individual from each related pair (the standard approach
    for unrelated GWAS).
-   Or model relatedness explicitly using mixed models if the study
    design includes family data.

Without this step, your GWAS results could be statistically biased and
less reproducible

#### 3. Population Genetics

Pi-hat helps investigate **population structure**, historical mixing,
and inbreeding patterns.

By summarizing how much relatedness exists within and between subgroups:

-   We can detect **hidden population clusters** or **cryptic
    substructure**.
-   We can study how migration and gene flow shape a population’s
    genetic landscape.
-   High average pi-hat values in a small population can indicate
    inbreeding, which is relevant for conservation genetics and human
    isolated populations.

In population studies, pi-hat complements other statistics like
F-statistics, kinship coefficients, or runs of homozygosity (ROH).

#### 4. Identifying Duplicate Samples

Pi-hat is also a **quality control tool** for catching errors like
**duplicate samples** or **sample swaps**.

A pi-hat value near **1.0** indicates that two samples are genetically
identical — they share nearly 100% of their genome IBD.

-   This may happen because the same DNA was genotyped twice under
    different IDs.
-   Or due to accidental swaps or contamination.

Detecting duplicates is critical for removing redundant or mislabeled
samples that would otherwise distort allele frequency estimates and
association results.

-   Check the relatedness. Use the independent SNPs (pruning) for this
    analysis and limit to autosomal chromosome only

#### Important Considerations When Using Pi-Hat Values

While pi-hat is a powerful and standard measure for estimating pairwise
genetic relatedness, it’s important to understand its **limitations**
and **context-specific challenges** to avoid misinterpretation.

#### 1. Inbreeding

In populations with **high levels of inbreeding**, individuals share a
larger fraction of their genome **IBD** with everyone in the population
— not just close relatives.

#### How does this affect pi-hat?

Pi-hat can be **inflated** for unrelated or distantly related pairs
because they all carry more shared segments from common ancestors. This
can **blur the line** between true close relatives and distant relatives
in an inbred population. For example, in small isolated communities or
animal breeding populations, unrelated individuals may appear more
related than they actually are.

#### What’s the risk?

If we don’t adjust for background inbreeding, we might:

-   Overestimate the number of close relatives.

-   Remove too many samples unnecessarily.

-   Misinterpret the degree of relatedness.

#### What to do:

Use additional measures, like F-statistics or runs of homozygosity
(ROH), to quantify inbreeding explicitly and adjust relatedness
estimates accordingly.

#### 2. Admixture

**Admixture** means individuals have mixed ancestry from **two or more
genetically distinct populations** (e.g., someone with European and
African ancestry).

#### How does this affect pi-hat?

Pi-hat assumes all individuals come from **a single homogeneous
population** with shared allele frequencies. But if individuals come
from mixed backgrounds, **allele frequency differences** between
populations can:

-   Cause pi-hat to **underestimate true relatedness** (e.g., a true
    sibling pair with highly admixed parents).
-   Or create **spurious relatedness** between individuals who just
    share similar ancestry proportions but aren’t actually closely
    related.

#### What’s the risk?

Admixture can **confound relatedness estimation**, leading to:

-   Missed detection of actual relatives.
-   False positives for related pairs.
-   Misleading inferences in population structure studies.

#### What to do:

Combine pi-hat with ancestry estimates (e.g., PCA, ADMIXTURE) to
interpret relationships in the context of known ancestry. Use tools that
can handle admixed data or adjust allele frequency estimates for
ancestry-specific components.

#### 3. Data Quality

The **accuracy** of pi-hat values depends directly on the **quality of
the genotype data**. Low-quality data, systematic genotyping errors, or
missing genotypes can distort IBD estimates.

#### How does this affect pi-hat?

-   **High missingness** → reduces the number of SNPs available for
    pairwise comparison → makes pi-hat less precise.
-   **Batch effects or poorly called SNPs** → may inflate or deflate
    shared segments.
-   **Non-random missingness** → can bias the estimates if some regions
    are missing in a correlated way across individuals.

#### What’s the risk?

Inaccurate pi-hat values can:

-   Miss real relatives.
-   Create false duplicate signals.
-   Lead to incorrect removal or retention of samples in QC.

#### What to do:

Always run robust **per-sample and per-SNP QC** first: filter
poor-quality SNPs, remove samples with high missingness, and check for
systematic batch effects. Use a high-quality set of common, well-called
SNPs when estimating relatedness.

### Pi-Hat Calculation

**Pi-hat *π̂*** estimates the **proportion of the genome** two
individuals share **IBD** — meaning they inherited segments from a
**recent common ancestor**.

When we compare two genomes, for every SNP we check:

-   Do they both have the same allele from the same ancestor (IBD)?

-   Or are their alleles identical by **state** (IBS) by coincidence?

The calculation estimates how much of the genome is **IBD**, adjusting
for the expected chance match if there was no relationship.

#### The Standard Formula

The core idea is:

$$
\large \hat\pi = P(IBD = 2) + \frac{1}{2} P(IBD=1)
$$

-   *P*(*I**B**D* = 2): Probability both alleles at a locus are IBD
    (e.g., identical twins).

-   *P*(*I**B**D* = 1): Probability one allele is IBD (e.g.,
    parent–child, full siblings).

-   *P*(*I**B**D* = 0): No alleles IBD.

So **pi-hat** sums the probability that:

-   Two individuals are identical at both alleles plus
-   Half the probability they share exactly one allele IBD.

#### Example Values

<table>
<thead>
<tr>
<th>Relationship</th>
<th>Expected pi-hat</th>
</tr>
</thead>
<tbody>
<tr>
<td>Identical twins or duplicate</td>
<td>1.0</td>
</tr>
<tr>
<td>Parent–offspring</td>
<td>0.5</td>
</tr>
<tr>
<td>Full siblings</td>
<td>0.5</td>
</tr>
<tr>
<td>Half siblings</td>
<td>0.25</td>
</tr>
<tr>
<td>First cousins</td>
<td>0.125</td>
</tr>
<tr>
<td>Unrelated</td>
<td>~0.0</td>
</tr>
</tbody>
</table>

#### How It’s Estimated in Practice

Tools like PLINK estimate pi-hat using genotype data:

1.  For each SNP, it looks at both individuals’ genotypes.
2.  It counts:

-   The proportion of SNPs where they share **0, 1, or 2** alleles IBS
    (identical by state).

1.  Using allele frequencies, it converts IBS counts into IBD estimates
    with a method-of-moments estimator.

The **method-of-moments** means it calculates the observed sharing and
adjusts for the expected sharing if individuals were unrelated.

#### PLINK output:

When we run `--genome`, PLINK calculates:

-   **Z0**: Estimated probability of IBD = 0
-   **Z1**: Estimated probability of IBD = 1
-   **Z2**: Estimated probability of IBD = 2

Then:

$$
\large \hat\pi = Z2 + \frac{1}{2}Z1
$$

#### PLINK command to keep only autosomes and exclude X chromosome

    .\plink     --bfile     2_QC_Raw_GWAS_data      --chr 1-22  --make-bed --out Autosomal

#### Autosomal QC checks assume diploidy

The above command is used to “Extract only autosomal SNPs (1–22) from my
raw GWAS data and save them as a new binary file called Autosomal.”

All of these previous steps, exclduing sex check, assume that every
chromosome is diploid. X chromosome breaks this rule.

#### It would distort heterozygosity calculations

-   Males have no heterozygosity on X chromosome SNPs in the non-PAR
    regions → this would artificially lower or skew overall
    heterozygosity if included.

-   So to check for **inbreeding or contamination**, we only use
    **autosomes**, where both sexes are comparable.

#### X chromosome needs separate handling

We do include the X chromosome later for:

-   Sex checks: Using X chromosome heterozygosity to verify reported
    sex.
-   Association testing: But with special handling (e.g., separate
    male/female coding).

So, in standard GWAS, **X chromosome QC is run in parallel but
separately** from the main autosomal pipeline.

-   Create independent SNPs through pruning

#### PLINK command to prune LD SNPs

    .\plink     --bfile     Autosomal   --indep-pairwise    50 5 0.2 --out raw-GWAS-data

-   \`.: This runs PLINK

-   `--bfile Autosomal`: Use my binary PLINK dataset named Autosomal

-   `--indep-pairwise 50 5 0.2`: This is the **LD pruning command** —
    one of **the most important GWAS QC steps** for creating **a set of
    independent SNPs** for downstream tasks like:

-   **PCA (principal component analysis)** to detect ancestry outliers.

-   **Relatedness** estimation.

-   Building a pruned set for kinship or mixed models.

    -   **50**: Use a sliding **window of 50 SNPs** at a time.
    -   **5**: Shift the **window by 5 SNPs** for each step.
    -   **0.2**: Within each window, if **any two SNPs have LD r² &gt;
        0.2, remove one to break the correlation**.

So this removes SNPs that are **highly correlated**, leaving you with a
**pruned set of approximately independent markers**.

-   This step will generate raw-GWAS-data.prune.in file.

#### Purpose:

LD pruning reduces redundancy due to linkage disequilibrium (LD). This
is important because:

-   **PCA**: LD can inflate signals and make ancestry PCs biased.
-   **IBD/Relatedness**: Using SNPs in high LD overestimates sharing.
-   **Mixed models**: LD-pruned sets make kinship estimates more
    accurate.

This step will generate a `raw-GWAS-data.prune.in` file

#### PLINK command

    .\plink --bfile 2_QC_Raw_GWAS_data -–extract raw-GWAS-data.prune.in  --genome --out related_check

-   `-–extract raw-GWAS-data.prune.in`: “Only use the pruned SNPs listed
    in raw-GWAS-data.prune.in.”
-   `--genome`: This tells PLINK to calculate pairwise IBD statistics
    for every possible pair of individuals in the dataset.

This step will generate `related_check.genome` to look for pairs with
pi-hat &gt; ~0.185. Decide which individual to keep from each pair
(usually the one with better QC, more complete data, or fewer missing
phenotypes).

    relatedness <- fread("D:/UNIX/GWAS/plink_linux_x86_64_20230116/Sex_check/Missing_Heter/Relatedness/related_check.genome")

    related2 <- subset(relatedness, relatedness$PI_HAT >=0.2)
    ggplot(relatedness, aes(x = PI_HAT, y = FID1))+
      geom_point(alpha=0.5, col = "#00bfff")+
      labs(x = "PI HAT value", y = "FID 1")+
      theme_classic()+
      scale_x_continuous(breaks = round(seq(0, 1, by = 0.1), 1))+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())+
      geom_vline(xintercept = 0.2, col = "Grey")+
      geom_point(data=relatedness %>%
                   filter(PI_HAT >= 0.2),
                 pch = 19,
                 size=1.6,
                 colour = "#e75480")+
      annotate(geom="text", x=0.27, y=1.0, label="Independent Samples",
               color="#003300")+
      annotate(geom="text", x=0.93, y=1.0, label="Duplicate",
               color="#003300")

<img src="Related_samples.png" alt="Relatedness" width="480" />
<p class="caption">
Relatedness
</p>

### Subsetting failed samples and creating a text file

    Failed_relatedness <- related2 %>% 
      select(1:2)

    ## Saving file 
    write.table(Failed_relatedness,  "Missing_Heter/Relatedness/relatedness_failed.txt", 
                row.names = FALSE, 
                col.names = FALSE,
                quote = FALSE)

Here `relatedness_failed.txt` contains failed related sample

#### Removing Failed samples from relatedness

    .\plink --bfile 2_QC_Raw_GWAS_data --remove relatedness_failed.txt --make-bed --out 3_QC_Raw_GWAS_data 

This step will create new binary files excluding failed related samples.

#### References

-   Purcell S et al., PLINK: A Tool Set for Whole-Genome Association and
    Population-Based Linkage Analyses, Am J Hum Genet, 2007.
-   Thompson EA (2013). Identity by descent: variation in meiosis,
    across genomes, and in populations. Genetics.
-   Weir BS, Cockerham CC (1984). Estimating F-Statistics for the
    Analysis of Population Structure. Evolution.
