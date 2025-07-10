# Identification of Individuals with elevated missing data rates

**Individual missingness** (sometimes called **sample call rate**) is
the proportion of genotypes that are **missing** for each person in the
dataset. In any GWAS, each individual is genotyped for hundreds of
thousands or millions of SNPs — but sometimes, due to technical or
biological reasons, not every SNP is successfully called for every
person.

### Why It Matters

High missingness in an individual’s data often signals **poor DNA
quality**, lab errors, or technical problems during genotyping (like low
DNA concentration, contamination, or plate failures). These low-quality
samples can introduce **noise and bias** in downstream association tests
— for example, by inflating false positives or causing subtle population
structure artifacts.

### What’s the Aim of This Check?

The goal is to **remove individuals** with too much missing data so
that: \* The remaining dataset is high-quality and complete. \* The
allele frequencies and test statistics are reliable. \* The chance of
spurious signals from poor-quality DNA is minimized.

#### Typical Threshold

A common threshold is to remove any sample with **&gt;2–5% missing
genotypes**. &gt; Very strict studies may use even tighter cutoffs, like
1%, especially for large datasets where small differences matter.

#### Practical Context

This check is always done **first** (or near the start) in QC pipelines,
because keeping bad samples distorts downstream SNP-level filters too —
e.g., if one person is missing 20% of their data, it can make some
perfectly fine SNPs look unreliable.

#### PLINK command to calculate missing rate.

    ./plink --bfile Sex_check_File -missing --out missing_data_rate

-   Command creates the files “missing\_data\_rate.imiss” and
    “missing\_data\_rate.lmiss”.
-   The “imiss” file (individual missingness) reports the proportion of
    missing genotypes per individual in the dataset. It lists each
    individual ID and the number and proportion of missing genotypes for
    that individual. This can be useful for identifying samples with
    high levels of missing data that may need to be removed from
    downstream analyses.

<table>
<thead>
<tr>
<th>FID</th>
<th>IID</th>
<th>MISS_PHENO</th>
<th>N_MISS</th>
<th>N_GENO</th>
<th>F_MISS</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>201320</td>
<td>N</td>
<td>631</td>
<td>536323</td>
<td>0.00177</td>
</tr>
<tr>
<td>1</td>
<td>201327</td>
<td>N</td>
<td>566</td>
<td>536323</td>
<td>0.00105</td>
</tr>
<tr>
<td>1</td>
<td>201419</td>
<td>N</td>
<td>94784</td>
<td>536323</td>
<td>0.1767</td>
</tr>
<tr>
<td>1</td>
<td>201567</td>
<td>N</td>
<td>16860</td>
<td>538448</td>
<td>0.03131</td>
</tr>
<tr>
<td>1</td>
<td>201359</td>
<td>N</td>
<td>1103</td>
<td>538448</td>
<td>0.00204</td>
</tr>
</tbody>
</table>

-   The fourth column in the .imiss file (N\_MISS) denotes the number of
    missing SNPs and the sixth column (F\_MISS) denotes the proportion
    of missing SNPs per individual.

Removing individuals with high missingness is one of the simplest but
most important QC steps in GWAS. It ensures that the dataset is built
only on reliable samples, protecting the validity of every association
result.

#### References

-   Anderson CA, Pettersson FH, Clarke GM, Cardon LR, Morris AP,
    Zondervan KT. (2010). Data quality control in genetic case-control
    association studies. Nat Protoc 5(9):1564–1573.
-   Marees AT, de Kluiver H, Stringer S, et al. (2018). A tutorial on
    conducting genome-wide association studies: Quality control and
    statistical analysis. Int J Methods Psychiatr Res. 27(2):e1608.
-   Purcell S, Neale B, Todd-Brown K, et al. (2007). PLINK: A Tool Set
    for Whole-Genome Association and Population-Based Linkage Analyses.
    Am J Hum Genet. 81(3):559–575. — Original paper describing PLINK’s
    implementation of missingness checks.
