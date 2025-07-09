A genome-wide association study (GWAS) is a type of study that uses
genomic data to identify genetic variations that are associated with a
particular trait or disease (Marees et al, 2018).The National Institute
of Health (NIH) defines GWAS as a study of common genetic variations
across the entire genome designed to identify genetic associations with
observable traits. A GWAS uses high-throughput genotyping technology. It
examines thousands of polymorphisms to relate them to a clinical
condition or measurable trait.It is a hypothesis-free approach to scan
markers across the whole genome.

It involves genotyping a large number of single nucleotide polymorphisms
(SNPs) across the entire genome of many individuals and comparing the
frequency of these SNPs in individuals with the trait or disease of
interest compared to those without the trait or disease.A traditional
**case-control study design** is the most common approach in GWASs. The
goal of a GWAS is to identify specific genetic loci that are associated
with the trait or disease, which can provide insight into the underlying
biological mechanisms and potentially lead to the development of new
diagnostic tools and treatments.

GWAS has been used to identify genetic risk factors for a wide range of
complex diseases and traits, including cardiovascular disease, type 2
diabetes, cancer, and many others. By utilizing large sample sizes and
high-throughput genotyping technologies, GWAS has revolutionized our
understanding of the genetic basis of complex traits and has become a
powerful tool in the field of genetics and genomics.

# Material and Methods

**Tools, Software and Files Required**

-   **plink**
    -   [link1](https://www.cog-genomics.org/plink2/)
    -   [link2](https://www.cog-genomics.org/plink/2.0/)
-   **bcftools and samtools**
    -   [link](http://www.htslib.org/download/)
-   **vcftools**
    -   [link](https://sourceforge.net/projects/vcftools/files/)
-   **Check for compatibility with the input required by the Sanger
    server.**
    -   [link](http://qbrc.swmed.edu/zhanxw/software/checkVCF/checkVCF-20140116.tar.gz)
-   **Haplotype Reference Consortium v1.1 panel for HRC site list**
    -   [link](HRC.r1-1.GRCh37.wgs.mac5.sites.tab)[link](ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/)
-   **HRC preparation checking tool**
    -   [link](https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9.zipa)
-   Set the path for these tools so we can directly call them from a
    different folder. To set the path first we have to enter in that
    specific directory

<!-- -->

    export PATH=$PATH:/mnt/d/UNIX/GWAS/plink_linux_x86_64_20230116
