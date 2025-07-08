# Per Individual Quality Control

**Per-sample QC** means checking the **quality of each individual’s
data** before running the genome-wide analysis. Its goal is to detect
and remove **bad samples** that could **introduce errors, bias, or false
associations**.

### Key Checks in Per-Sample QC

-   **Missingness**: Samples with a high proportion of missing genotype
    calls are likely to have poor DNA quality or technical errors and
    should be removed.

-   **Sex Check**: Confirms that the genetically inferred sex (from X
    chromosome data) matches the reported sex, helping catch sample
    swaps or labelling mistakes.

-   **Heterozygosity**: Samples with unusually high or low
    heterozygosity may indicate contamination (mixed DNA) or inbreeding
    — both can distort results.

-   **Relatedness**: GWAS assumes that individuals are unrelated.
    Closely related individuals inflate association statistics, so
    duplicates or relatives must be identified and usually removed.

-   **Ancestry Check**: Principal component analysis (PCA) is used to
    detect samples whose genetic ancestry does not match the target
    population, which could cause spurious associations due to
    population structure.

If poor-quality samples are not removed, they can lead to **false
signals**, inflate test statistics, or hide true genetic effects. Robust
per-sample QC makes sure your study is **reliable, replicable, and
scientifically valid**.

Though there are many steps in QC. But one of the main step could be
converting ped and map file into binary files.

### Converting ped and map into binary format

#### PLINK command

    ./plink –file raw_GWAS_data --make-bed --out plink

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
<td style="text-align: left;"><code>--out plink</code></td>
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

# References

1- Marees, A.T., et al, 2018. A tutorial on conducting genome‐wide
association studies: Quality control and statistical analysis. *Int J
Methods Psychiatr Res*, Jun; 27(2): e1608.

2- Anderson, C.A. et al, 2010. Data quality control in genetic
case-control association studies. *Nat Protoc*, Sep:5(9):1564-73

3- Singh, Sandeep Kumar, “A Case-Only Genome-wide Association Study of
Gender- and Age-specific Risk Markers for Childhood Leukemia” (2015).
FIU Electronic Theses and Dissertations. 1832
