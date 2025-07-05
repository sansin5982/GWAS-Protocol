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
<td>0.027</td>
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
    -   PEDSEX
    -   SNPSEX
    -   STATUS (OK, PROBLEM)
    -   F

-   The above command also provide a log file, here
    GWAS\_Sex\_Check.log, that provides many information including
    number of cases and controls, males and females count, individuals
    with ambiguous code, etc

-   Extract the IDs of individuals with discordant sex information. In
    situations where discrepancies cannot be resolved, remove the
    individuals through following command. <br> <br>

**PLINK command to remove the individuals based on sex information**
<br> <br> **plink –bfile raw\_GWAS\_data –remove
discordant-sex-individuals-file.txt –make-bed –out
1\_QC\_Raw\_GWAS\_data** <br> <br> (File
“discordant-sex-individuals-file.txt”, should contain only FID and IID
of the individuals that have to be removed)

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
