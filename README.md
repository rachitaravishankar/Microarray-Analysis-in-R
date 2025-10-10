# Microarray Analysis in R

Analysis of the ImmGen Microarray Phase 2 dataset (GSE37448) from Mus musculus, exploring gene-expression profiles of immune cells using Affymetrix Mouse Gene 1.0 ST arrays.

# Dataset Overview

- Title: ImmGen Microarray Phase 2
- Organism: Mus musculus
- Platform: Affymetrix Mouse Gene 1.0 ST Array (GPL6246)
- Samples: 194
- Experiment type: Expression profiling by array
- Status: Public (April 24, 2012)
- Citation:
Elpek KG, Cremasco V, Shen H, Harvey CJ et al. The tumor microenvironment shapes lineage, transcriptional, and functional diversity of infiltrating myeloid cells.
Cancer Immunol Res. 2014 Jul;2(7):655–67. PMID: 24801837

# Summary

This dataset was generated as part of the Immunological Genome Project (ImmGen).
It provides gene-expression microarray data using Ambion WT Expression Kit reagents instead of the Affymetrix GeneChip WT cDNA Synthesis and Amplification Kits.
The data has been used to explore how the tumor microenvironment shapes transcriptional and functional diversity in infiltrating myeloid cells.

# Data Access

GSE37448 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37448

# Planned Analysis
Analysis is performed on two datasets : series matrix file and dataset deposited in the Immgen website.

For the raw series matrix file the following steps are followed - 
- Quality control of CEL files
- Normalization (RMA)
- Differential expression analysis (limma)
- Enrichment Analysis (Gene Ontology based)

For the immgen dataset the following steps are followed - 
- Dimensionalty Reduction

# References

ImmGen Project Website

NCBI GEO: GSE37448

Elpek KG, Cremasco V, Shen H, et al. The tumor microenvironment shapes lineage, transcriptional, and functional diversity of infiltrating myeloid cells. Cancer Immunol Res. 2014;2(7):655-667. doi:10.1158/2326-6066.CIR-13-0209

*Disclaimer: This project is intended solely for personal use and practice. It is not an official analysis or reproduction of the ImmGen project.*
