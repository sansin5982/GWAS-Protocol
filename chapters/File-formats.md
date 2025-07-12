# File format to perform PLINK

## Basic input files

We expect two basic input files, Raw.ped and Raw.map, for association
analysis. Here, raw is the file name. PED (Plink Pedigree File) and MAP
(Plink Map File) are file formats commonly used in genome-wide
association studies (GWAS). These files can also be converted into
binary format (.bed, .bim, .fam).

-   Covariate file for variables like age, sex, or principal components.

#### Initial Checks:

-   Verify file integrity and consistency (e.g., matching sample IDs
    across genotype, phenotype, and covariate files).
-   Ensure all datasets use the same genome build (e.g., GRCh 37 or 38)
    to avoid positional mismatches.

### PED File

The PED file is a tab-delimited text file that contains information
about individuals in a study, including their family relationships and
genotype data. Each row in the PED file represents a single individual,
and each column represents a different type of information, such as the
individual’s family ID, individual ID, and genotype data. The first six
columns in a PED file are required and include the following
information:

<table>
<thead>
<tr>
<th>FID</th>
<th>IID</th>
<th>PID</th>
<th>MID</th>
<th>Sex</th>
<th>P</th>
<th>rs1</th>
<th>rs2</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>1</td>
<td>CT</td>
<td>AG</td>
</tr>
<tr>
<td>2</td>
<td>2</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>CG</td>
<td>AA</td>
</tr>
<tr>
<td>3</td>
<td>3</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>0</td>
<td>CC</td>
<td>TT</td>
</tr>
</tbody>
</table>

-   FID = Family ID,a unique identifier for the family to which the
    individual belongs
-   IID = Individual ID, a unique identifier for the individual within
    the family
-   PID = Paternal ID,the individual ID of the father (-9 if unknown)
-   MID = Maternal ID, the individual ID of the mother (-9 if unknown)
-   Sex (1 = male; 2 = female; other = -9)
-   P = Phenotype, the phenotype of the individual, usually coded as a
    binary trait (e.g. affected=1, unaffected=2, missing=0)

The remaining columns in the PED file contain genotype data for each
individual.

-   rs1 = SNP1
-   rs2 = SNP2

Sixth column represents phenotype. The phenotype can be either an
affection status or quantitative trait. Plink can automatically detect
the phenotype based on the code.

### MAP file

Each line of the MAP file explains a single marker and contains four
columns. The MAP file is a tab-delimited text file that contains
information about each SNP (single nucleotide polymorphism) that was
genotyped in the study. Each row in the MAP file represents a single
SNP, and each column represents a different type of information, such as
the chromosome number, SNP identifier, and base-pair position. The
columns in a MAP file include the following information:

<table>
<thead>
<tr>
<th>Chr</th>
<th>SNP</th>
<th>GD</th>
<th>BPP</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>rs1</td>
<td>0</td>
<td>870050</td>
</tr>
<tr>
<td>1</td>
<td>rs2</td>
<td>0</td>
<td>870150</td>
</tr>
<tr>
<td>1</td>
<td>rs3</td>
<td>0</td>
<td>870322</td>
</tr>
<tr>
<td>1</td>
<td>rs4</td>
<td>0</td>
<td>870878</td>
</tr>
</tbody>
</table>

-   Chr = Chromosome (1-22, X, Y or 0 if unplaced), the chromosome
    number on which the SNP is located
-   rs = a unique identifier for the SNP
-   GD = the distance between this SNP and the previous SNP (in Morgans)
-   BPP = the position of the SNP on the chromosome (in base-pairs)

The autosomes should be coded 1 through 22. The following other codes
can be used to specify other chromosome types: <br>

<table>
<thead>
<tr>
<th>Region</th>
<th>Genomic region</th>
<th>Number</th>
</tr>
</thead>
<tbody>
<tr>
<td>X</td>
<td>X chromosome</td>
<td>23</td>
</tr>
<tr>
<td>Y</td>
<td>Y chromosome</td>
<td>24</td>
</tr>
<tr>
<td>XY</td>
<td>Pseudo-autosomal region of X</td>
<td>25</td>
</tr>
<tr>
<td>MT</td>
<td>Mitochondrial</td>
<td>26</td>
</tr>
</tbody>
</table>

Together, the PED and MAP files provide all the information needed to
perform a GWAS. The genotype data in the PED file can be combined with
the information about each SNP in the MAP file to test for associations
between genotype and phenotype. This is typically done using a
statistical software package, such as PLINK.

## Binary PED files

Binary file formats are an alternative to the text-based PED/MAP file
formats commonly used in GWAS. Binary formats offer several advantages
over text-based file formats, including faster loading times, reduced
storage space requirements, and improved data compression. Binary files
contain extended map file (a.bim) having information about allele names,
binary ped file (a\*bed) and pedigree/phenotype information in separate
file (a.fam).

### BED

The BED file contains information on the presence or absence of a
specific allele for each marker at each sample, represented as binary
values (0 or 1).

### BIM

A BIM file is a file format used to store and organize data on genetic
variants that are included in a study. The BIM file contains information
on the chromosome, variant ID, position, and allele information for each
genetic variant.

The BIM file is typically a tab-delimited text file with six columns:
chromosome number, variant ID, genetic position, base-pair position, the
reference allele, and the alternate allele. The chromosome number column
specifies which chromosome the genetic variant is located on, while the
variant ID column provides a unique identifier for the variant. The
genetic position and base-pair position columns provide information on
the location of the variant on the chromosome. The reference allele and
alternate allele columns provide information on the two possible alleles
for the variant.

GWAS typically involve the analysis of large numbers of genetic
variants, and the BIM file helps to ensure that data is organized and
analyzed correctly. The BIM file is used in conjunction with other
files, such as genotype files (e.g., PLINK binary files), to identify
associations between genetic variants and traits or diseases of
interest.

<table>
<thead>
<tr>
<th>Chr</th>
<th>SNP</th>
<th>GD</th>
<th>BPP</th>
<th>Allele1</th>
<th>Allele2</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>rs1</td>
<td>0</td>
<td>870050</td>
<td>C</td>
<td>T</td>
</tr>
<tr>
<td>1</td>
<td>rs2</td>
<td>0</td>
<td>870150</td>
<td>A</td>
<td>G</td>
</tr>
<tr>
<td>1</td>
<td>rs3</td>
<td>0</td>
<td>870322</td>
<td>G</td>
<td>G</td>
</tr>
<tr>
<td>1</td>
<td>rs4</td>
<td>0</td>
<td>870878</td>
<td>C</td>
<td>A</td>
</tr>
</tbody>
</table>

### FAM

In a GWAS, a FAM file is a text file that contains information about the
individuals being studied. The FAM file is used to specify the familial
relationships and phenotype information for each individual in the
study.

The FAM file typically contains three columns of data separated by
whitespace or tab characters. The first column contains a family ID,
which is a unique identifier for a particular family or group of related
individuals. The second column contains an individual ID, which is a
unique identifier for each individual within a family. The third column
contains the phenotype information for each individual, such as disease
status or a quantitative trait.

For example, a FAM file might look like this:

<table>
<thead>
<tr>
<th>familyID</th>
<th>individualID</th>
<th>fatherID</th>
<th>motherID</th>
<th>sex</th>
<th>phenotype</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>1</td>
</tr>
<tr>
<td>1</td>
<td>2</td>
<td>0</td>
<td>0</td>
<td>2</td>
<td>0</td>
</tr>
<tr>
<td>1</td>
<td>3</td>
<td>1</td>
<td>2</td>
<td>1</td>
<td>1</td>
</tr>
<tr>
<td>1</td>
<td>4</td>
<td>1</td>
<td>2</td>
<td>2</td>
<td>0</td>
</tr>
</tbody>
</table>

In this example, there is one family with four individuals. The first
column specifies the family ID, which is 1 in this case. The second
column specifies the individual ID, ranging from 1 to 4. The third and
fourth columns specify the IDs of the father and mother for each
individual, respectively. The fifth column specifies the sex of each
individual (1 for male and 2 for female). The sixth column specifies the
phenotype, with 1 indicating affected by the disease or trait of
interest, and 0 indicating unaffected.The phenotype can be either a
quantitative trait or an affection status column: PLINK will
automatically detect which type.

## What is a VCF File?

### VCF = Variant Call Format

A **VCF** file is a **standard text format** developed by the **1000
Genomes Project** to store **genetic variant data** (like SNPs and
indels) along with metadata and genotypes.

It’s widely used in:

-   GWAS (Genome-Wide Association Studies)
-   Whole-genome/exome sequencing
-   Genotype imputation
-   Reference panel sharing (like 1000 Genomes or HRC)

### Structure of a VCF File

There are two main parts:

#### 1. Header lines (start with \##)

Contain meta-information, such as:

    ##fileformat=VCFv4.2
    ##source=GATK
    ##contig=<ID=1,length=248956422>

#### 2. Column header line (starts with \#CHROM)

    #CHROM  POS  ID     REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2 ...

#### 3. Data rows

Example:

    1  752566  rs3094315  C  T  .  PASS  .  GT  0/0  0/1

This line represents:

-   Chromosome 1
-   Position 752566
-   rsID = rs3094315
-   REF allele = C
-   ALT allele = T
-   Genotype: Sample1 is homozygous reference (0/0), Sample2 is
    heterozygous (0/1)

### Why VCF is Important for GWAS

-   Stores **all variants** (SNPs, indels, multi-allelic)
-   Encodes **sample genotypes** for thousands of individuals
-   Used for **data exchange** (e.g., imputation servers, reference
    panels)
-   Flexible and extensible

VCFs are used **before GWAS**, especially after sequencing or
imputation, to be converted into binary PLINK format for analysis.

### Conversion: VCF to PLINK Format

Tools Required:

-   PLINK
-   bcftools
-   vcftools

We want to convert this:

    mydata.vcf.gz

To PLINK’s binary format:

    mydata.bed, mydata.bim, mydata.fam

Step-by-step:

    # Make sure your VCF is bgzipped and indexed
    bgzip mydata.vcf
    tabix -p vcf mydata.vcf.gz

    # Convert to PLINK binary format
    plink --vcf mydata.vcf.gz --make-bed --out mydata

This will produce:

-   `mydata.bed`: genotype matrix
-   `mydata.bim`: SNP info
-   `mydata.fam`: sample info

### Convert PLINK to VCF

Use this command:

    plink --bfile mydata --recode vcf --out mydata

It will generate:

-   `mydata.vcf` (uncompressed)

We can compress/index it afterward:

    bgzip mydata.vcf
    tabix -p vcf mydata.vcf.gz

#### Considerations During Conversion

<table>
<colgroup>
<col style="width: 18%" />
<col style="width: 81%" />
</colgroup>
<thead>
<tr>
<th>Issue</th>
<th>Consideration</th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>Genome build</strong></td>
<td>Must match the build used in reference panel or imputation</td>
</tr>
<tr>
<td><strong>rsID missing</strong></td>
<td>Make sure your <code>.bim</code> file has correct rsIDs before
converting</td>
</tr>
<tr>
<td><strong>Multiallelics</strong></td>
<td>VCF handles these; PLINK expects biallelic — use
<code>--biallelic-only strict</code> if needed</td>
</tr>
<tr>
<td><strong>Sex chromosomes</strong></td>
<td>PLINK requires careful handling of X/Y during conversion</td>
</tr>
</tbody>
</table>

### Combining Multiple Cohorts

    plink --bfile dataset1 --bmerge dataset2.bed dataset2.bim dataset2.fam --make-bed --out merged

-   Resolve strand issues or SNP mismatches using `--flip` or
    `--exclude` if merge errors occur.

#### References

-   PLINK: <https://www.cog-genomics.org/plink/>
-   VCF: <https://samtools.github.io/hts-specs/VCFv4.3.pdf>
