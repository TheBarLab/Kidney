# Kidney

Scripts used to generate the data for "Kidney-specific methylation patterns correlate with kidney function and are lost upon kidney disease progression" by Naor Sagy, Noa Meyrom, Pazit Beckerman, Oren Pleniceanu and Daniel Z Bar.

20220914_unique_cg_kidney.csv - list of kidney unique methylation sites and their values across various tissues. Based on NGDC-CNCB.

Confounder_Adjustment_code.py - script used to calculate confounder effect on R and p-values (Supp. Table 2). 

CorrelationDist_IF.py - script used to calculate correlation distribution (Supp. Fig.  2). 

FindUniqueSites.py - script used to identify tissue-unique methylation sites for all tissues.

HM450.hg19_min_data.xlsx - sites coordinates and chromosomes, based on the Illumina manifesto.

IF-eGFR.py - script for IF-eGFR correlation (Supp. Fig. 8).

Methyl_sites_heatmap_V1.R - R script for creating figure 5 heatmap.

PullOutSequences.py - script for pulling DNA sequences from Entrez.

chinese_v1_metadata.csv - metadata for the NGDC-CNCB dataset used.

kidney_unique_sites_0_2.csv - output example for unique sites.
