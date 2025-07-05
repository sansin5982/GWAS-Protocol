<script type="text/javascript" async
    src="https://polyfill.io/v3/polyfill.min.js?features=es6">
</script>
<script type="text/javascript" async
    src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.0/es5/tex-mml-chtml.js">
</script>

# Per Individual Quality Control

Per-individual QC screens genotype to identify subjects that may
introduce bias, if not removed. Poorly genotyped (low call rate)
individuals will increase error in the study and may significantly
affect the results. There are several steps of per-individual QC for a
GWAS data set.

### Step 1: Converting ped and map into binary format

#### PLINK command

    ./plink –file raw_GWAS_data --make-bed --out raw_GWAS_data

#### What each part means

<table>
<colgroup>
<col style="width: 47%" />
<col style="width: 52%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Part</th>
<th>Meaning</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><code>./plink</code></td>
<td>Runs the PLINK executable from the current directory
(<code>./</code>).</td>
</tr>
<tr>
<td style="text-align: left;"><code>--file raw_GWAS_data</code></td>
<td>Tells PLINK to use <strong>text format PED/MAP files</strong> named
<code>raw_GWAS_data.ped</code> and <code>raw_GWAS_data.map</code>.</td>
</tr>
<tr>
<td style="text-align: left;"><code>--make-bed</code></td>
<td>Converts the input files into <strong>binary PLINK format</strong>:
<code>.bed</code> / <code>.bim</code> / <code>.fam</code>.</td>
</tr>
<tr>
<td style="text-align: left;"><code>--out raw_GWAS_data</code></td>
<td>To give a new name to file</td>
</tr>
</tbody>
</table>

#### what happens here?

-   The `--file` option means: &gt; “Use the old text-based PED/MAP
    format.”

-   The `--make-bed` option means: &gt; “Create binary files instead.”

So this step **converts** our dataset from **PED/MAP** → **BED/BIM/FAM
format**.

#### Why do this?

The **binary format** is:

-   Faster to read/write
-   Smaller file size
-   Required for many downstream PLINK operations

Most people do **all GWAS work** in binary format (`--bfile`).

#### What files are created?

PLINK will produce:

-   plink.bed → genotype binary data
-   plink.bim → variant info
-   plink.fam → sample info

By default they are named `plink.*` unless we add `--out`:

### Step 2: Identification of individuals with discordant sex information

This was the first step of QC, performed to identify subjects that have
inconclusive/contradictory gender information, as it can lead to
spurious associations. In a case-control study, samples that show the
wrong gender information are suggested to be excluded from further QC
and statistical analysis. Subjects were appropriately recoded or
removed, if information was inconclusive, for further analyses.

#### PLINK command

    ./plink --bfile raw_GWAS_data --check-sex --out GWAS_Sex_Check

#### What each part means

<table>
<colgroup>
<col style="width: 55%" />
<col style="width: 44%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Part</th>
<th style="text-align: left;">Meaning</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><code>./plink</code></td>
<td style="text-align: left;">Run the PLINK software executable from the
current directory (<code>./</code> means “this folder”).</td>
</tr>
<tr>
<td style="text-align: left;"><code>--bfile raw_GWAS_data</code></td>
<td style="text-align: left;">Use the <strong>binary PLINK
files</strong>: <code>raw_GWAS_data.bed</code>,
<code>raw_GWAS_data.bim</code>, and <code>raw_GWAS_data.fam</code>.
These 3 files together define your genotype dataset.</td>
</tr>
<tr>
<td style="text-align: left;"><code>--check-sex</code></td>
<td style="text-align: left;">Run PLINK’s <strong>sex check</strong>
routine. This checks whether the genetically inferred sex matches the
sex reported in your <code>.fam</code> file.</td>
</tr>
<tr>
<td style="text-align: left;"><code>--out GWAS_Sex_Check</code></td>
<td style="text-align: left;">Name the output files: e.g.,
<code>GWAS_Sex_Check.sexcheck</code> will be created.</td>
</tr>
</tbody>
</table>

#### What does –check-sex actually do?

-   PLINK uses **X chromosome data** to estimate the inbreeding
    coefficient **F** for each sample:

    -   If F ≈ 1 → likely male (only one X chromosome → more
        homozygosity)
    -   If F ≈ 0 → likely female (two X chromosomes → more
        heterozygosity)

-   It then compares this estimate with the **sex** listed in your
    **.fam** file (column 5).

-   If there’s a mismatch, it flags it:

-   Possible reasons: sample swap, data entry error, or sex chromosome
    abnormality.

#### Output

-   The main output is:
    -   GWAS\_Sex\_Check.sexcheck

<table>
<thead>
<tr>
<th>FID</th>
<th>IID</th>
<th>PEDSEX</th>
<th>SNPSEX</th>
<th>STATUS</th>
<th>F</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>201320</td>
<td>2</td>
<td>2</td>
<td>OK</td>
<td>0.017</td>
</tr>
<tr>
<td>1</td>
<td>201320</td>
<td>1</td>
<td>1</td>
<td>OK</td>
<td>0.976</td>
</tr>
<tr>
<td>1</td>
<td>201327</td>
<td>2</td>
<td>1</td>
<td>PROBLEM</td>
<td>0.974</td>
</tr>
<tr>
<td>1</td>
<td>201335</td>
<td>2</td>
<td>0</td>
<td>PROBLEM</td>
<td>0.456</td>
</tr>
<tr>
<td>1</td>
<td>201342</td>
<td>1</td>
<td>2</td>
<td>PROBLEM</td>
<td>0.632</td>
</tr>
<tr>
<td>1</td>
<td>201359</td>
<td>1</td>
<td>0</td>
<td>PROBLEM</td>
<td>0.321</td>
</tr>
</tbody>
</table>

-   Each row has:

    -   FID
    -   IID
    -   PEDSEX: The reported sex in our `.fam` file
    -   SNPSEX: The inferred sex estimated from the genotype data
    -   STATUS (OK, PROBLEM)
    -   F: Give value. If it is &lt; 0.2, sample is suggested as female
        and if score is &gt; 0.8 sample is represented as male.

-   The above command also provide a log file, here
    `GWAS_Sex_Check.log`, that provides many information including
    number of cases and controls, males and females count, individuals
    with ambiguous code, etc

-   Extract the IDs of individuals with discordant sex information. In
    situations where discrepancies cannot be resolved, remove the
    individuals through following command.

#### PLINK command to remove the individuals based on sex information

    plink --bfile raw_GWAS_data --remove discordant-sex-individuals-file.txt --make-bed --out 1_QC_Raw_GWAS_data

File `“discordant-sex-individuals-file.txt”`, should contain only FID
and IID of the individuals that have to be removed)

<table>
<thead>
<tr>
<th>FID</th>
<th>IID</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>201320</td>
</tr>
<tr>
<td>1</td>
<td>201327</td>
</tr>
<tr>
<td>1</td>
<td>201335</td>
</tr>
<tr>
<td>1</td>
<td>201342</td>
</tr>
<tr>
<td>1</td>
<td>201359</td>
</tr>
</tbody>
</table>

    Gender <- read.table("Sex_check_1.sexcheck", header = T, as.is = T) %>%
      na.omit()
    png("Gender_check.png")

    ggplot(Gender, aes(x=F, y= PEDSEX, col = STATUS))+
      geom_point()+
      labs(y="Gender", x = "F score")+
      theme_classic()+
      scale_x_continuous(breaks = round(seq(min(Gender$F), 1, by = 0.01), 1))+
      theme(strip.text.x = element_blank())+
      geom_circle(aes(x0 = 0.61,  y0 = 2, r = 0.42),
                  inherit.aes = FALSE,
                  col = "Red")+
      geom_circle(aes(x0 = 0.05,  y0 = 1, r = 0.46),
                  inherit.aes = FALSE,
                  col = "Red")+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none")+
      annotate(geom="text", x=-0.2, y=2.1, label="Females (< 0.2)",
               color="#e75480")+
      annotate(geom="text", x=0.9, y=1.1, label="Males (> 0.8)",
               color="#00bfff")+
      annotate(geom="text", x=0.61, y=2.1, label="Females failed (> 0.2)",
               color="#e75480")+
      annotate(geom="text", x=0.05, y=1.1, label="Males failed (< 0.8)",
               color="#00bfff")
    dev.off()

The above R script will create a png file containing information based
on STATUS column and differentiate males, females and ambiguous data
separately.

<img src="Gender_check.png" alt="Discordant Sex information" width="480" />
<p class="caption">
Discordant Sex information
</p>

The above figure explains Individuals Discordant sex information.
Samples failed QC are in red circle. Some samples also showed value less
than 0. It is good to cross check these samples or **exclude it F score
is &lt; -0.05**.

### Step 3: Identification of individuals with elevated missing data rates

This QC step is used to identify and exclude individuals with too much
missing genotype data. Genotype accuracy and genotype call rate can be
significantly affected by variations in DNA quality. A high genotype
failure rate suggest poor DNA sample quality.

#### PLINK command to calculate missing rate.

    ./plink2 --bfile 1_QC_Raw_GWAS_data --missing --out missing_data_rate

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
<td>N1103</td>
<td>538448</td>
<td>0.00204</td>
<td></td>
</tr>
</tbody>
</table>

-   The fourth column in the .imiss file (N\_MISS) denotes the number of
    missing SNPs and the sixth column (F\_MISS) denotes the proportion
    of missing SNPs per individual.

### Step 4: Identification of individuals with outlying heterozygosity rate

#### PLINK command

    ./plink --bfile 1_QC_Raw_GWAS_data --het --out outlying_heterozygosity_rate

<table>
<thead>
<tr>
<th>FID</th>
<th>IID</th>
<th>O(HOM)</th>
<th>E(HOM)</th>
<th>N(NM)</th>
<th>F</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>201320</td>
<td>228640</td>
<td>2.259e+05</td>
<td>325867</td>
<td>0.0272</td>
</tr>
<tr>
<td>1</td>
<td>201337</td>
<td>227550</td>
<td>2.26e+05</td>
<td>326067</td>
<td>0.01506</td>
</tr>
</tbody>
</table>

**NOTE: plink2 format will give results in a different way**

-   Command creates the file “outlying\_heterozygosity\_rate.het”, in
    which the third column denotes the observed number of homozygous
    genotypes \[O(Hom)\] and the fifth column denotes the number of
    non-missing genotypes \[N(NM)\] per individual.

<!-- -->

    # Missing individual & Heterozygosity rate
    miss <- fread("SEX_data/Missing_sample/missing_data_rate.imiss")
    hetro <- fread("SEX_data/Missing_sample/Heterozygosity_rate/outlying_heterozygosity_rate.het")
    head(miss, 2)
    head(hetro, 2)

Next we calculate observed heterozygosity rate

$$
\text {Observed Heterozygosity Rate} = \frac {\text {non-missing genotypes (N(NM)) - Observed number of homozygous genotypes (O(HOM))}}{\text {non-missing genotypes (N(NM))}} 
$$

    # Calculate the observed heterozyosity rate
    hetro$obs_hetero_rate <- ((hetro$`N(NM)`)-hetro$`O(HOM)`)/hetro$`N(NM)`

-   Merge the miss and hetro dataframe above created

<!-- -->

    # Merge missing file and heterozygoisty file
    hetro_miss <- miss %>% 
      left_join(hetro, by = "IID")

    # Creating plot
    png("Missing_hetero_check.png")
    ggplot(hetro_miss, aes(x = F_MISS, y = obs_hetero_rate))+
      geom_point(alpha = 0.5, col = "#00bfff")+
      labs(x ="Proportion of Missing genotypes(log scale)", y = "Heterozygosity rate")+
      scale_x_log10(limits = c(0.0001, 1))+
      theme_classic()+
      scale_y_continuous(limits = c((min(hetro_miss$obs_hetero_rate) - 0.02), (max(hetro_miss$obs_hetero_rate) + 0.02)))+
      geom_hline(yintercept = mean(hetro_miss$obs_hetero_rate) + 3*sd(hetro_miss$obs_hetero_rate), col = "Grey")+
      geom_hline(yintercept = mean(hetro_miss$obs_hetero_rate) - 3*sd(hetro_miss$obs_hetero_rate), col = "Grey")+
      geom_vline(xintercept = 0.006, col = "Grey")+
      geom_point(data=hetro_miss %>%
                   filter(F_MISS > 0.01),
                 pch = 19,
                 size=1.6,
                 colour = "#e75480")+
      annotate(geom="text", x=0.035, y=0.303, label="Missing genotypes (>1%)",
               color="#003300")+
      annotate(geom="text", x=0.1, y=0.285, label="Excess heterozygosity rate (Â± 3 sd from mean)",
               color="#003300")+
      annotate(geom="text", x=0.00017, y=0.315, label="+3 sd from mean",
               color="#003300")+
      annotate(geom="text", x=0.00017, y=0.301, label="-3 sd from mean",
               color="#003300")+
      annotate(geom="text", x=0.0015, y=0.31, label="Genotyping > 99%",
               color="#003300")

    dev.off()

The above R script will create a png file containing information based
on missing rate and heterozygosity rate. For missing rate we set a
threshold value &lt; 0.01 and for heterzygosity rate between + and - 3
standard deviation.

<img src="Missing_hetero_check.png" alt="Individual missingness and heterozygoisty rate" width="480" />
<p class="caption">
Individual missingness and heterozygoisty rate
</p>

### Step 5: Identification of duplicate samples

In a population based study, it is important that all samples should be
unrelated. Presence of duplicate, first- or second-degree relatives will
introduce bias in the study as their genotypes will be overrepresented.
This step was used to identify all related and duplicate individuals for
removal. A metric (identity by state, IBS) for each pair of individuals
was calculated to identify duplicate samples. IBS is defined as, at a
locus, two individuals who have an identical nucleotide sequence or the
same allele.

<table>
<thead>
<tr>
<th>PIHAT values</th>
<th>Relation</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>Identical or Duplicate</td>
</tr>
<tr>
<td>0.8 and above</td>
<td>Highly related</td>
</tr>
<tr>
<td>0.5</td>
<td>First Degree</td>
</tr>
<tr>
<td>0.25</td>
<td>Highly related</td>
</tr>
<tr>
<td>0.125</td>
<td>First Cousin</td>
</tr>
<tr>
<td>0.0625</td>
<td>Second Cousin</td>
</tr>
</tbody>
</table>

-   Check the relatedness. Use the independent SNPs (pruning) for this
    analysis and limit to autosomal chromosome only

#### PLINK command

    .\plink     --bfile     2_QC_Raw_GWAS_data      --chr 1-22  --make-bed --out Autosomal

-   Create independent SNPs through pruning

#### PLINK command

    .\plink     --bfile     Autosomal   --indep-pairwise    50 5 0.2 --out raw-GWAS-data

it will generate raw-GWAS-data.prune.in file. This file use in next step
\* Check relatedness

#### PLINK command

    plink --bfile 2_QC_Raw_GWAS_data –extract raw-GWAS-data.prune.in  --genome --out related_check

<img src="Related_samples.png" alt="Relatedness" width="480" />
<p class="caption">
Relatedness
</p>
