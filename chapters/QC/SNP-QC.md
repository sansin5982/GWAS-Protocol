## Per SNP quality control

#### Identification of SNPs with elevated missing data rates

#### Test markers for different genotype call rates between cases and contols

#### Minor Allele Frequency

#### SNPs failing Hardy Weinberg equilibrium

-   If the frequency of observed genotypes of a variant in a population
    can be derived from the observed allele frequencies, the genetic
    variant is said to be in Hardy–Weinberg equilibrium.
-   Genotype frequencies remain stable from one generation to another in
    the absence of any evolutionary pressure (selection, mutation and
    migration).
-   The goodness-of-fit test is used to test HWE is not a powerful test
    and reliability of results depends on the sample size.

<table>
<thead>
<tr>
<th>Sample size</th>
<th>50</th>
<th>100</th>
<th>200</th>
<th>300</th>
<th>400</th>
</tr>
</thead>
<tbody>
<tr>
<td>AA count</td>
<td>21</td>
<td>42</td>
<td>84</td>
<td>126</td>
<td>168</td>
</tr>
<tr>
<td>AB count</td>
<td>25</td>
<td>50</td>
<td>100</td>
<td>150</td>
<td>200</td>
</tr>
<tr>
<td>BB count</td>
<td>4</td>
<td>8</td>
<td>16</td>
<td>24</td>
<td>32</td>
</tr>
<tr>
<td>HWE p value</td>
<td>0.52</td>
<td>0.26</td>
<td>0.08</td>
<td><strong>0.02</strong></td>
<td><strong>0.009</strong></td>
</tr>
</tbody>
</table>

Source: Genetic Epidemiology: Mehmet T Dorak

-   Rule of thumb: the heterozygote frequency can only reach a maximum
    of 50%. If heterozygote frequencies are more than 50%, it is a clear
    sign of HWD, regardless of statistical test result.
-   Most common reason is not biological, **genotyping error** is most
    plausible exploration.

#### PLINK Command

    ./plink2 --bfile 4_QC_Raw_GWAS_data --geno 0.01 --hwe 0.00000001 --make-bed --out 5_QC_Raw_GWAS_data

#### 
