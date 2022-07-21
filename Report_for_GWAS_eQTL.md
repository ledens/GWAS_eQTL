# GWAS_eQTL Report
R code for co-localization of GWAS and eQTL candidates for the purpose of selecting and prioritizing genes and SNPs

Irritable bowel syndrome (IBS) is a dysregulation of the gut immune system that affects an estimated 11% of the global population. (1)  IBS has a multifactorial pathophysiology that includes changes to the intestinal microbiome and intestinal permeability, along with mucosal inflammation and immune system activation. (1) Twin studies suggest autoimmune diseases have a roughly 50% genetic basis, with environmental factors being responsible for the other half. Peters et al. (2) suggests that analysis of gene expression across different immune cell types is important for finding associations between genotype and phenotype. (2) Aschenbrenner et al. (3) identifies monocyte cells as one of the sources of interleukin (IL)-23, a cytokine that plays a central role in the inflammatory bowel response. (3) In this report, I perform an association between IBS genome-wide association study (GWAS) data and monocyte cell expression quantitative trait loci (eQTL) data for the purpose of selecting and prioritizing genes and single nucleotide polymorphisms (SNPs).

## Approach 

### GWAS

GWAS data for the phenotype “irritable bowel syndrome” was obtained under the relevant trait page (EFO_0000555) of the GWAS Catalog. (4) This catalog used the genome assembly GRCh38.p13 and dbSNP Build 154. All association studies were downloaded to a single database and explored using R. The data consisted of 28 variant entries. 15 of the entries came from a single study, Bonfiglio et al. (5), while a further 7 came from a second meta-analysis study by the same authors. (6) The primary Bonfiglio et al. (5) study focused on SNP associations on chromosome 9 for 9,576 IBS patients (336,499 controls). All the patients in this study were women from Europe and the United States between the ages of 40 and 69 years. As this study is the primary contributor to my GWAS data, my results will be biased toward this population type.

Of the 28 GWAS SNPs, 12 were intergenic and 16 were intragenic. For the intragenic SNPs, 1 was identified as a missense variant, 2 were non-coding exon variants, and the rest were intron variants. This means that an active role played by the GWAS IBS variants will be focused around gene regulation. 

### eQTLs

eQTL data for monocyte cells was downloaded from the Database of Immune Cell Expression (DICE). (8) Specifically, the eQTLs were downloaded in vcf format from the filtered Monocyte, classical database (Build date 6.15.18, DICE Build ID: DICE-DB 1). This data was filtered to only contain eQTLs with P-values less than 1 x 10-4 and consisted of 315,456 entries. As this cutoff was a lower P-value than those mentioned by Doerge and Churchill (9) and Farris et al. (10), no additional P-value filter was applied. All entries were given as single-location SNPs without coordinates ranges.

Problematically, the DICE catalog used the genome assembly GRCh37.p19. To translate genome assemblies, I initially attempted to use the R package ‘rsnps’ (11), which provides in-R access to data from OpenSNP and NCBI’s dbSNP SNP databases. This would have allowed me to convert all SNP coordinates to the same genome build. However, the online database call within rsnps is time-gated and in trial runs took roughly 1-3 seconds per SNP to return a result. As this method was too slow, I instead employed the online Ensembl Assembly Converter (12) to convert the genome coordinates of my entire eQTL dataset to genome assembly GRCh38. Using this online converter first required a slight modification to the DICE file format, which was accomplished with the command sed -i 's/chr//' . Unfortunately, Assembly Converter only converted the file’s chromosome coordinates, and therefore the gene associations within the eQTL data were no longer usable. 

### Co-localization

Co-location between GWAS and eQTL data was determined by chromosome coordinate proximity. This method exclusively focused on cis-eQTL associations and did not locate any trans-eQTLs. The GWAS SNPs were divided into an intergenic and an intragenic set. A different co-localization strategy was used for each set.

Intragenic GWAS SNPs were associated with a mapped gene. Manual checks with the Ensembl Genome Browser (13) ensured the GWAS mapped genes agreed with the chromosome coordinates. All genes that contained at least one GWAS SNP were selected. There was one SNP that occurred in two genes, and therefore both genes were considered. This method assumes any SNP marker within the gene is highlighting an IBS association of that entire gene. Therefore, the final rankings focused on the identified genes and not on the individual SNPs. Gene coordinate ranges were found by using the biomaRt Ensembl annotation Bioconductor R package. (14) eQTL SNP locations were then checked to determine how many occurred within each gene. Final gene rankings were determined by a combination of internal eQTL count, eQTL average P-value, the identifying SNP P-value, and the identifying SNP risk allele value.  The following genes were found to be high-priority: STK17A, DAGLB and EPSTI1. STK17A, also known as DRAK (13), was identified to have the highest priority, and is involved in death-associated protein kinase activity. (15) 

One problem with the above intragenic GWAS SNPs method is it limits the eQTL search to occurrences inside the gene. In this instance, no eQTL SNP peaks were observed interior to any selected gene. However, if the window is expanded to +/- the average gene size (~0.3 Mb) on either side, then peaks can be seen nearby for the genes ABLIM2, EPSTI1, MTSS1, and STK17A. If eQTL peak association is selected as the sole ranking factor, then ABLIM2 becomes the highest priority gene. ABLIM2 produces an actin binding LIM protein, which is part of a group of proteins that play in several biological processes, such as cancer differentiation. (16)

Intergenic GWAS SNPs were not associated with any genes. Instead, a +/- basepair window centered on the SNP location was used to check if any eQTL SNPs were located inside the window. The window size was determined relative to the average basepair range of the genes identified by the intragenic GWAS SNPs. Several window sizes were sampled, with the size increasing until every intergenic GWAS SNP had at least one associated eQTL. The intergenic GWAS SNP and eQTL association was explored in detail for a full window size of ~0.6 Mb. While 7 GWAS SNPs were associated with eQTLs at this window size, most contained only a few eQTL associations. 

Only SNPs rs2523599 and rs4663866 are near any significant numbers of eQTLs. Interestingly, while rs2523599 has more eQTLs SNPs within its window, rs4663866 is clearly closer to an eQTL peak. However, the eQTL peak has a low P-value, so I would still prioritize rs2523599 over rs4663866. 

## Caveats

•	53% of my data came from a Europe and United States women-only study. This will bias my final rankings toward that population group.

•	Accepting low P-value GWAS data greatly increased my false discovery odds.

•	My eQTL data is point-based SNPs, without estimated ranges, meaning I am relying only on the SNP markers when determining associations based on location.

•	I focused only on one cell type (monocyte cells) due to time constraints (more cell types were originally planned)

•	Using a coordinate-based approach means I only looked for cis-eQTLs.

## Discussion

There are multiple ways to perform the co-localization of GWAS and eQTL candidates to aid in prioritization. I applied two related approaches to look at co-localization between IBS GWAS data and monocyte cell eQTLs. The first approach used intragenic GWAS SNPS to identify genes of interest. Then eQTLs inside those genes were identified. The interior eQTL count and the GWAS SNP P-value and risk allele frequency were used in combination to prioritize the selected genes. Alternatively, if eQTL identification was expanded to chromosomal regions nearby the gene, better eQTL co-localization could be achieved. If this improved co-localization was considered as the primary prioritization factor, then the gene rankings changed. This alternative prioritization better matches the second co-localization approach, but does not include factors such as the GWAS P-value or risk allele frequency and is therefore inferior to the previous ranking system. The second co-localization approach used intergenic GWAS SNPS to identify chromosomal regions of interest. Then eQTLs inside those regions were identified after an optimal region size was determined. Only two intergenic GWAS SNPS were found to be adjacent to eQTL peaks and were prioritized based on their proximity.

Overall, the largest issues I faced with my co-localization came from my GWAS data. Not only is the available data limited, but it is derived primarily from a single source. Additionally, my data had large relative P-values, meaning my false discovery rates are likely also very large. If a co-localization between IBS GWAS and eQTL candidates was to be repeated, I would recommend more exploratory methods that could potentially identify GWAs susceptible regions. 

## References

1.	Fairbrass KM, Costantino SJ, Gracie DJ, Ford AC. Prevalence of irritable bowel syndrome-type symptoms in patients with inflammatory bowel disease in remission: a systematic review and meta-analysis. The Lancet Gastroenterology & Hepatology. 2020;5(12):1053-62.
2.	Peters JE, Lyons PA, Lee JC, Richard AC, Fortune MD, Newcombe PJ, et al. Insight into Genotype-Phenotype Associations through eQTL Mapping in Multiple Cell Types in Health and Immune-Mediated Disease. PLOS Genetics. 2016;12(3):e1005908.
3.	Aschenbrenner D, Quaranta M, Banerjee S, Ilott N, Jansen J, Steere B, et al. Deconvolution of monocyte responses in inflammatory bowel disease reveals an IL-1 cytokine network that regulates IL-23 in genetic and acquired IL-10 resistance. Gut. 2021;70(6):1023-36.
4.	The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019 [Cited Dec 13, 2021]. 2019. Available from: https://www.ebi.ac.uk/gwas/efotraits/EFO_0000555.
5.	Bonfiglio F, Zheng T, Garcia-Etxebarria K, Hadizadeh F, Bujanda L, Bresso F, et al. Female-Specific Association Between Variants on Chromosome 9 and Self-Reported Diagnosis of Irritable Bowel Syndrome. Gastroenterology. 2018;155(1):168-79.
6.	Bonfiglio F, Henström M, Nag A, Hadizadeh F, Zheng T, Cenit MC, et al. A GWAS meta-analysis from 5 population-based cohorts implicates ion channel genes in the pathogenesis of irritable bowel syndrome. Neurogastroenterol Motil. 2018;30(9):e13358.
7.	Chen Z, Boehnke M, Wen X, Mukherjee B. Revisiting the genome-wide significance threshold for common variant GWAS. G3 Genes|Genomes|Genetics. 2021;11(2).
8.	Ha B, Greenbaum JA, Shmiedel BJ, Singh D, Madrigal A, Valdovino-Gonzalez AG, et al. Database of Immune Cell EQTLs, Expression, Epigenomics. The Journal of Immunology. 2019;202(1 Supplement):131.18-.18.
9.	Doerge RW, Churchill GA. Permutation Tests for Multiple Loci Affecting a Quantitative Character. Genetics. 1996;142(1):285-94.
10.	Farris SP, Wolen AR, Miles MF. Using expression genetics to study the neurobiology of ethanol and alcoholism. Int Rev Neurobiol. 2010;91:95-128.
11.	The rsnps package. rOpenSci. 2021 [cited Dec. 6, 2021]. Available from: https://docs.ropensci.org/rsnps/.
12.	Assembly converter Ensembl: European Bioinformatics Institute - European Molecular Biology Laboratory; 2021 [Cited Dec 13, 2021]. Available from: https://uswest.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core.
13.	Howe KL, Achuthan P, Allen J, Allen J, Alvarez-Jarreta J, Amode MR, et al. Ensembl 2021. Nucleic Acids Research. 2020;49(D1):D884-D91.
14.	Smith ML, Durinck S, Huber W. Accessing Ensembl annotation with biomaRt Bioconductor. 2021 [Cited Dec 13, 2021]. 2.50.1. Available from: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
15.	Sanjo H, Kawai T, Akira S. DRAKs, novel serine/threonine kinases related to death-associated protein kinase that trigger apoptosis. J Biol Chem. 1998;273(44):29066-71.
16.	Krupp M, Weinmann A, Galle PR, Teufel A. Actin binding LIM protein 3 (abLIM3). Int J Mol Med. 2006;17(1):129-33


