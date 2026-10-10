---
title: "Chapter 3 — Genetic Variation and Genetic Architecture"
output:
  md_document:
    variant: markdown
    preserve_yaml: true
---

## Learning objectives

By the end of this chapter, we should be able to distinguish major
classes of DNA variation; interpret reference and alternate alleles
without confusing them with risk alleles; explain allele frequency and
minor-allele frequency; calculate genotype and allele counts by hand;
distinguish effect size from statistical significance; describe
monogenic, oligogenic, and polygenic architectures; and recognize how
linkage disequilibrium and study design affect the interpretation of
GWAS signals.

**Scope.** We focus on the biological and statistical meaning of genetic
variation. Population-genetic mechanisms, Hardy--Weinberg equilibrium,
detailed linkage disequilibrium, and population stratification are
developed in Chapter 4. Genotyping technology and file formats are
covered in later chapters. All numerical examples in this chapter are
**hypothetical teaching examples**, not results from our Hapmap dataset.

## 3.1 What is genetic variation?

Human DNA sequences are highly similar, but they are not identical. A
**genetic variant** is a difference in DNA sequence relative to a stated
reference sequence or between individuals. A **locus** is a genomic
position or region; an **allele** is one of the sequence states observed
at that locus. The term *variant* can refer to a single-base change, an
insertion or deletion, or a much larger rearrangement.

We should distinguish three questions:

1.  **What changed?** The variant class and sequence representation.
2.  **How common is it?** The frequency of the allele in a defined
    population or sample.
3.  **What does it do?** Its functional or phenotypic effect, which may
    be unknown.

A variant can be common without being biologically important, or rare
with a large effect. Frequency alone does not establish pathogenicity.

## 3.2 Major classes of genetic variation

### Single-nucleotide variants and SNPs

A **single-nucleotide variant (SNV)** differs at one nucleotide
position, for example A versus G. A **single-nucleotide polymorphism
(SNP)** is commonly used for a single-base variant present in a
population; historically, *polymorphism* often implied a frequency of at
least 1%, although contemporary usage is not completely uniform. When
describing a specific call without assuming population frequency,
**SNV** is the more general term.

At an autosomal biallelic SNP with alleles A and G, the diploid
genotypes can be AA, AG, and GG. The genotype is a pair of alleles; it
is not the same thing as an allele frequency.

### Insertions and deletions

An **insertion** adds sequence relative to the reference; a **deletion**
removes sequence. Small insertions and deletions are called **indels**.
For example, `ACGT` versus `AGT` differs by a deleted C in one
representation. In repetitive DNA, the same biological event can be
written at more than one coordinate until variants are normalized
against a specific reference genome.

### Structural variation

**Structural variants (SVs)** generally involve larger changes, often
operationally defined as at least 50 base pairs, although conventions
vary. They include deletions, duplications, inversions, insertions, and
translocations. **Copy-number variants (CNVs)** alter the number of
copies of a DNA segment. Not all structural variants change copy number:
an inversion can be copy-neutral.

### Repeated sequences and complex variation

**Short tandem repeats (STRs)** consist of repeated short sequence
motifs, such as `CAG` repeated multiple times. Repeat expansions can
affect disease biology. Some genomic regions contain complex or
multiallelic variation that cannot be adequately described as one
biallelic SNP. Classical HLA alleles are an important example of highly
polymorphic, sequence-defined alleles; SNPs in the MHC region are not
interchangeable with classical HLA allele calls.

::: {.figure style="text-align: center"}
`<img src="figures/ch03-fig-3-1-1.png" alt="Figure 3.1. Schematic examples of sequence-variant classes. Sequences are invented for teaching and are not mapped to a genomic coordinate."  />`{=html}
<p class="caption">
Figure 3.1. Schematic examples of sequence-variant classes. Sequences
are invented for teaching and are not mapped to a genomic coordinate.
</p>
:::

**Important limitation:** The figure illustrates sequence differences;
it does not provide VCF-compliant, left-normalized variant
representations.

## 3.3 Reference, alternate, effect, and minor alleles

Several allele labels answer different questions and should never be
treated as synonyms.

  -----------------------------------------------------------------------
  Term                    Meaning                 Key caution
  ----------------------- ----------------------- -----------------------
  Reference allele (REF)  Sequence matching the   Not necessarily the
                          chosen reference        common or ancestral
                          assembly at the         allele
                          represented locus       

  Alternate allele (ALT)  Sequence differing from Not necessarily the
                          REF                     risk allele

  Effect allele           Allele to which a       Its choice determines
                          reported effect         the sign of a
                          estimate refers         regression coefficient

  Other allele            Comparator for the      May be REF or ALT
                          effect allele in a      
                          biallelic analysis      

  Minor allele            Less frequent allele in Can change between
                          a specified             populations
                          population/sample       

  Ancestral allele        Inferred allele state   Requires evidence; not
                          in an evolutionary      synonymous with REF
                          ancestor                
  -----------------------------------------------------------------------

For example, if a GWAS reports an effect allele G with an odds ratio of
1.20 relative to A, the association is reported **per additional G
allele** under the stated model. This does not tell us whether G is REF,
ALT, ancestral, or minor.

### Why allele harmonization matters

Suppose one study reports the effect of G versus A, and another reports
A versus G at the same biallelic SNP. Under consistent model
assumptions, the effect coefficient changes sign when the effect allele
is switched:

$$\beta_{A}=-\beta_{G}.$$

The corresponding odds ratio becomes its reciprocal:

$$OR_{A}=\frac{1}{OR_{G}}.$$

If $OR_{G}=1.20$, then $OR_{A}=1/1.20\approx0.833$. This is a change of
**reference direction**, not a disagreement between studies. Strand
alignment and ambiguous A/T or C/G SNPs require additional care in real
data.

## 3.4 Counting alleles and genotypes

Consider 100 unrelated diploid individuals typed at a biallelic
autosomal SNP:

  Genotype      Number of individuals   Number of G alleles contributed
  ----------- ----------------------- ---------------------------------
  AA                               36                                 0
  AG                               48                                48
  GG                               16                                32
  **Total**                   **100**                            **80**

Each individual contributes two allele copies, so the total number of
observed allele copies is $2(100)=200$. The frequency of G is:

$$ f_G=\frac{n_{AG}+2n_{GG}}{2N}
     =\frac{48+2(16)}{2(100)}
     =\frac{80}{200}=0.40.$$

The frequency of A is $f_A=1-f_G=0.60$. Because G is less frequent in
this particular sample, its **minor-allele frequency (MAF)** is 0.40.
Its **minor-allele count (MAC)** is 80. The genotype frequencies are
0.36, 0.48, and 0.16, respectively.

For a **biallelic** site, $MAF=\min(f_A,f_G)$. This shortcut should not
be applied without modification to multiallelic variants. If genotypes
are missing, the denominator must reflect the number of successfully
observed allele copies, not automatically twice the full study sample
size.

``` r
genotype_counts <- c(AA=36, AG=48, GG=16)
N_called <- sum(genotype_counts)
G_count <- genotype_counts["AG"] + 2 * genotype_counts["GG"]
G_frequency <- unname(G_count / (2 * N_called))
A_frequency <- 1 - G_frequency
MAF <- min(G_frequency, A_frequency)
MAC <- min(G_count, 2*N_called - G_count)
data.frame(N_called=N_called, G_count=unname(G_count),
           G_frequency=G_frequency, A_frequency=A_frequency,
           MAF=MAF, MAC=unname(MAC))
```

    ##   N_called G_count G_frequency A_frequency MAF MAC
    ## 1      100      80         0.4         0.6 0.4  80

``` r
stopifnot(isTRUE(all.equal(G_frequency, 0.4)))
```

### Frequency is population-specific

An allele frequency describes a particular sample or population at a
particular time. A variant with frequency 0.02 in one population may be
substantially more frequent in another. Thus, labels such as *rare* or
*common* must be interpreted with the reference population and frequency
threshold in mind.

A widely used descriptive convention is **common** (MAF at least 5%),
**low-frequency** (roughly 1% to below 5%), and **rare** (below 1%).
These are conventions, not biological laws, and some studies use
different cutoffs. For a sample of 100 diploid individuals, an allele
present once has MAC 1 and an observed frequency of $1/200=0.005$
(0.5%).

::: {.figure style="text-align: center"}
`<img src="figures/ch03-fig-3-2-1.png" alt="Figure 3.2. Hypothetical genotype counts and the corresponding allele copies at one biallelic SNP."  />`{=html}
<p class="caption">
Figure 3.2. Hypothetical genotype counts and the corresponding allele
copies at one biallelic SNP.
</p>
:::

## 3.5 From genotype to an association effect

GWAS typically tests whether genotype is associated with a phenotype.
For an additive autosomal model, we can code the number of G alleles as
$G=0,1,2$ for AA, AG, and GG, respectively. This is called **allele
dosage**. The choice of G as effect allele is a modeling convention.

For a quantitative trait such as systolic blood pressure, a simple
illustrative linear model is:

$$Y=\beta_0+\beta_GG+\varepsilon,$$

where $Y$ is the trait, $\beta_0$ is the predicted value at $G=0$,
$\beta_G$ is the predicted change per additional G allele, and
$\varepsilon$ captures variation not explained by this simple model.

Suppose $\beta_0=120$ mmHg and $\beta_G=2$ mmHg. The predicted means
are:

  Genotype     G dosage   Predicted systolic blood pressure
  ---------- ---------- -----------------------------------
  AA                  0                 $120+2(0)=120$ mmHg
  AG                  1                 $120+2(1)=122$ mmHg
  GG                  2                 $120+2(2)=124$ mmHg

This example assumes an additive effect on the **measurement scale** and
omits covariates. It does not imply that every individual with GG has
exactly 124 mmHg or that the variant causes the observed difference.

::: {.figure style="text-align: center"}
`<img src="figures/ch03-fig-3-3-1.png" alt="Figure 3.3. Predicted trait means under a hypothetical additive model. These are model predictions, not observed participant measurements."  />`{=html}
<p class="caption">
Figure 3.3. Predicted trait means under a hypothetical additive model.
These are model predictions, not observed participant measurements.
</p>
:::

### Additive, dominant, and recessive models

An additive model uses dosage 0, 1, 2. A **dominant model** for G
compares AG/GG with AA (coding 0, 1, 1). A **recessive model** for G
compares GG with AA/AG (coding 0, 0, 1). These are statistical genotype
codings; the word *dominant* in a regression model should not be taken
as proof of a classical Mendelian inheritance pattern.

  Genotype     Additive G   Dominant G   Recessive G
  ---------- ------------ ------------ -------------
  AA                    0            0             0
  AG                    1            1             0
  GG                    2            1             1

``` r
genotypes <- c("AA", "AG", "GG")
data.frame(genotype=genotypes,
           additive=c(0,1,2), dominant=c(0,1,1), recessive=c(0,0,1))
```

    ##   genotype additive dominant recessive
    ## 1       AA        0        0         0
    ## 2       AG        1        1         0
    ## 3       GG        2        1         1

In large GWAS, additive models are often a practical starting point
because they use one per-variant dosage coefficient. More complex models
can be useful when biologically justified, but they require appropriate
sample sizes and interpretation.

## 3.6 Effect size, significance, and precision

An **effect size** describes the magnitude and direction of an estimated
association. For quantitative traits, it may be a regression coefficient
in trait units or standardized units. For binary traits, logistic
regression commonly reports an **odds ratio (OR)**. A **P value**
quantifies compatibility of the observed test statistic with a specified
null model; it is not the probability that a biological hypothesis is
true.

An OR of 1.10 represents 10% higher **odds** per effect allele under the
model, not necessarily a 10% increase in disease probability or risk. An
OR of 0.90 represents 10% lower odds. Confidence intervals summarize
uncertainty in the effect estimate, conditional on the model and
assumptions.

A tiny effect can achieve a small P value in a very large sample; a
substantial estimated effect can be imprecise in a small sample.
**Statistical significance and biological importance are different
questions.** Detailed inference and association-model mathematics will
be covered in Chapters 6 and 37--42.

## 3.7 What is genetic architecture?

The **genetic architecture** of a phenotype describes how genetic
variants collectively contribute to variation in that phenotype.
Important dimensions include the number of contributing loci,
effect-size distribution, allele frequencies, modes of action,
interactions, and the extent to which variants are correlated.

### Monogenic, oligogenic, and polygenic traits

A **monogenic disorder** is primarily attributable to pathogenic
variation in one gene, although penetrance, expressivity, and genetic
background may still vary. An **oligogenic** model involves a small
number of contributing genes or loci. A **polygenic** trait is
influenced by many genetic variants, often each with a small effect.

Many complex phenotypes---including height, blood pressure, and
susceptibility to common psychiatric disorders---have highly polygenic
components. However, describing a trait as polygenic does not mean that
every associated variant is causal or that environmental influences are
unimportant.

::: {.figure style="text-align: center"}
`<img src="figures/ch03-fig-3-4-1.png" alt="Figure 3.4. Conceptual genetic architectures. The number and sizes of plotted effects are illustrative, not estimates for a particular phenotype."  />`{=html}
<p class="caption">
Figure 3.4. Conceptual genetic architectures. The number and sizes of
plotted effects are illustrative, not estimates for a particular
phenotype.
</p>
:::

### Penetrance and expressivity

**Penetrance** is the proportion of individuals with a specified
genotype who exhibit a defined phenotype under stated conditions.
**Expressivity** describes variation in the phenotype's manifestations
among individuals with that genotype. These concepts are particularly
useful for Mendelian disorders but also remind us that genotype does not
mechanically determine every observed outcome.

### Heritability is not determinism

**Heritability** is a population- and environment-dependent measure of
how much variation in a phenotype is statistically attributable to
genetic differences under a particular model. It does not mean that a
given percentage of an individual's phenotype is genetically determined.
Nor does high heritability imply that environmental interventions cannot
alter the phenotype.

GWAS-derived SNP heritability usually captures only a specified
component of genetic variation tagged by the analyzed markers under the
estimation method. It is not automatically identical to family-based
heritability. Detailed estimation belongs in specialized post-GWAS
methods rather than this introductory chapter.

## 3.8 Why an associated SNP may not be causal

Nearby variants can be inherited together more often than expected by
chance, a phenomenon called **linkage disequilibrium (LD)**. A genotyped
SNP may therefore show association because it is correlated with another
variant that influences the phenotype. In such a case, the observed SNP
is a **tagging marker**, not necessarily the causal variant.

::: {.figure style="text-align: center"}
`<img src="figures/ch03-fig-3-5-1.png" alt="Figure 3.5. Conceptual distinction between a tested tagging SNP and a possible causal variant. The association is illustrative; the diagram does not establish causality."  />`{=html}
<p class="caption">
Figure 3.5. Conceptual distinction between a tested tagging SNP and a
possible causal variant. The association is illustrative; the diagram
does not establish causality.
</p>
:::

LD will be introduced quantitatively in Chapter 4. Later, fine-mapping
and functional studies can help prioritize candidate causal variants,
but those techniques are addressed in our separate advanced genomics
resources.

## 3.9 Coding, noncoding, and regulatory variants

A **coding variant** falls within a protein-coding sequence and may be
synonymous, missense, nonsense, or otherwise alter translation. A
**noncoding variant** lies outside protein-coding sequence; it may
affect transcription, splicing, chromatin regulation, RNA function, or
have no detectable effect. A variant's location near a gene does not
prove that gene is the effector.

Functional interpretation requires attention to genome build, transcript
choice, gene annotation, tissue context, and supporting evidence. GWAS
association by itself cannot tell us which transcript or molecular
mechanism is affected.

## 3.10 Practical R exercise: frequencies and dosage

We will simulate a small **teaching-only** set of genotypes. This is not
a simulation of real Hapmap samples, and no PLINK QC threshold is being
applied.

``` r
set.seed(303)
N <- 200
p_G <- 0.30
# Independent allele draws under a simplified illustrative model.
allele_1 <- rbinom(N, 1, p_G)
allele_2 <- rbinom(N, 1, p_G)
G_dosage <- allele_1 + allele_2
observed_G_frequency <- sum(G_dosage) / (2*N)
counts <- table(factor(G_dosage, levels=0:2,
                       labels=c("AA", "AG", "GG")))
print(counts)
```

    ## 
    ## AA AG GG 
    ## 94 86 20

``` r
cat("Observed G frequency:", round(observed_G_frequency, 4), "\n")
```

    ## Observed G frequency: 0.315

``` r
cat("Expected input frequency:", p_G, "\n")
```

    ## Expected input frequency: 0.3

``` r
stopifnot(sum(counts)==N, all(G_dosage %in% 0:2))
```

The observed frequency need not equal 0.30 exactly because the example
uses random draws. Our simulation assumes independent allele draws and
is not a realistic model of every population or genomic region. In
particular, it does not model LD, relatedness, population structure, or
genotyping error.

::: {.figure style="text-align: center"}
`<img src="figures/ch03-fig-3-6-1.png" alt="Figure 3.6. Genotype distribution from the reproducible simulated R example (N = 200; input G-allele frequency = 0.30)."  />`{=html}
<p class="caption">
Figure 3.6. Genotype distribution from the reproducible simulated R
example (N = 200; input G-allele frequency = 0.30).
</p>
:::

### Methods and reproducibility record

  -----------------------------------------------------------------------
  Item                                Record
  ----------------------------------- -----------------------------------
  Input                               No external genotype files; values
                                      generated within this chapter

  Software                            Base R; record the local R version
                                      when rendering

  Random seed                         `303`

  Simulated sample size               200 diploid individuals

  Assumed G-allele frequency          0.30

  Generated data                      Two independent Bernoulli allele
                                      draws per individual

  Verification                        All dosages in 0--2; genotype
                                      counts sum to 200

  QC exclusions                       None; this is not a QC analysis

  Limitations                         Does not represent real Hapmap
                                      observations, LD, ancestry, or
                                      errors
  -----------------------------------------------------------------------

``` r
R.version.string
```

    ## [1] "R version 4.5.3 (2026-03-11 ucrt)"

## 3.11 Common misconceptions

- **"ALT means risk allele."** Incorrect. ALT is defined relative to a
  reference assembly; risk direction comes from a phenotype-specific
  association model.
- **"The minor allele is always the same worldwide."** Incorrect. Allele
  frequencies differ between populations.
- **"A significant GWAS SNP is the causal mutation."** Not necessarily.
  The marker may tag a causal variant through LD or reflect confounding.
- **"A polygenic trait is entirely genetic."** Incorrect. Environmental
  and measurement factors may contribute substantially.
- **"A large odds ratio always means a large absolute risk
  difference."** Incorrect. Absolute risk also depends on baseline risk
  and model context.

## Chapter summary

Genetic variation ranges from single-base substitutions to large
structural rearrangements. Reference, alternate, effect, minor, and
ancestral alleles describe different properties. For a biallelic diploid
locus, we can calculate allele frequency from genotype counts and code
genotypes as 0, 1, or 2 effect-allele copies. Genetic architecture
concerns the number, frequency, and effects of contributing variants.
Many common traits are polygenic, and an associated marker is not
necessarily causal. These distinctions are essential before we study
population genetics and the statistical assumptions underlying GWAS.

## Concept checks and exercises

1.  Explain the difference between an SNV, an indel, and a copy-number
    variant.
2.  Why is an ALT allele not automatically the minor allele or the risk
    allele?
3.  In 150 successfully genotyped diploid individuals, we observe 75 AA,
    60 AG, and 15 GG. Calculate the G-allele count, G-allele frequency,
    and MAF. **Check:** G count = 90; frequency = 0.30; MAF = 0.30.
4.  If a study reports $OR_G=1.25$, what is the odds ratio when A
    becomes the effect allele? **Check:** $OR_A=0.80$.
5.  Write the additive, dominant-G, and recessive-G codes for AA, AG,
    and GG.
6.  Explain why a large GWAS sample may identify a variant with a small
    effect.
7.  What does it mean for a trait to be polygenic? Why does this not
    establish genetic determinism?
8.  Why can a tagging SNP show association without being causal?
9.  Explain why the sample minor allele can change when the study
    population changes.
10. Modify the R simulation to use `p_G <- 0.10`. Re-run with the same
    seed and compare the genotype distribution; label the results as
    simulated.

## References and further reading

1.  Visscher PM, Wray NR, Zhang Q, et al. 10 years of GWAS discovery:
    Biology, function, and translation. *American Journal of Human
    Genetics*. 2017;101:5--22.
    <https://doi.org/10.1016/j.ajhg.2017.06.005>
2.  Abdellaoui A, Yengo L, Verweij KJH, Visscher PM. 15 years of GWAS
    discovery: Realizing the promise. *American Journal of Human
    Genetics*. 2023;110:179--194.
    <https://doi.org/10.1016/j.ajhg.2022.12.011>
3.  Marees AT, de Kluiver H, Stringer S, et al. A tutorial on conducting
    genome-wide association studies: Quality control and statistical
    analysis. *International Journal of Methods in Psychiatric
    Research*. 2018;27:e1608. <https://doi.org/10.1002/mpr.1608>
4.  National Human Genome Research Institute. Talking Glossary of
    Genomic and Genetic Terms.
    <https://www.genome.gov/genetics-glossary>
5.  The 1000 Genomes Project Consortium. A global reference for human
    genetic variation. *Nature*. 2015;526:68--74.
    <https://doi.org/10.1038/nature15393>

**Next chapter:** Chapter 4 --- Population Genetics: allele-frequency
differences, Hardy--Weinberg equilibrium, genetic drift, linkage
disequilibrium, and the foundations of population structure.
