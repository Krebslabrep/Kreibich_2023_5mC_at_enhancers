# Kreibich et al., 2023 - The regulation of enhancers by DNA methylation

# Introduction
This GitHub repository contains all primary code to reproduce the main analyses for the publication "Single molecule footprinting identifies context-dependent regulation of enhancers by DNA methylation" [Kreibich et al., 2022, bioRxiv](https://doi.org/10.1101/2022.05.19.492653), from the [Krebs laboratory](https://www.embl.de/research/units/genome_biology/krebs/index.html).

In this publication, we use Single Molecule Footprinting (SMF) to measure chromatin accessibility (CA), transcription factor (TF) binding and endogenous DNA methylation (5mC) simultaneously at the same single DNA molecules. This enables us to analyze the molecular relationship between these gene regulatory features.  
SMF is a high-throughput sequencing technology developed in the [Krebs laboratory](https://www.embl.de/research/units/genome_biology/krebs/index.html). It consists of marking accessible genomic cytosines in the GpC context using a exogenous methytransferase and subsequent bisulfite sequencing (BS). Consequently, GpCs protected by the binding of DNA-interacting proteins (e.g. TFs, nucleosomes, GTFs, etc.) will remain unmethylated, while the accessible cytosines will be methylated.
In this project, we read out the cytosine methylation in CpG and GpC context independently for each DNA molecule at genomic bins of interest (101 bp bins or around TF binding sites (TFBS)). Those methylation scores are then used to determine the association between DNA methylation (CpG context) and chromatin accessibility or TF binding (GpC context) using a statistical test (Cochran-Mantel-Haenszel test).  

Other analysis scripts in this repository include the analysis of Precision Run-On data (PRO-seq) and ChIP-seq data.  

# SMF Analyses strategies
To exclude ambigious cytosine methylation calls (e.g. GCGs), we use more stringent cytosine contexts:  
NWCGW - for CpGs   
DGCHN - for GpCs

In general, the SMF analysis is performed with two strategies: (1) Chromatin accessibility view and (2) TFBS view.

**(1) Chromatin accessibility (CA) analysis**
* This analysis is based on the assumption that the binding of a TF to the DNA opens up a wider region that is marked as accessible.
* Here, we use a 101 bp bin surrounding 10-60% methylated CpG to calculate the local chromatin accessibility using the mean of the GpC methylation calls.
* We then perform a Cochran-Mantel-Haenszel test to associate the CpG methylation of the center CpG with the accessible or inaccessible fraction.

A more detailed schematic view of this analysis can be found in Figure S1H of the manuscript and is also shown here:  
![Kreibich et al, 2023. Figure S1H - Schematic description of CA analysis strategy](Figure_S1H_Analysis_scheme.png)

**(2) Transcription factor binding site (TFBS) analysis**
* Here, we center a 30 bp center bin at TFBS that contains at least one CpG. In addition 2 neighboring X bp bins up- and downstream are created. The GpC methylation within tose three bins is used to determine the binding state of the TFBS (TF bound, accessible, nucleosome bound).
* We then perform a Cochran-Mantel-Haenszel test to associate the CpG methylation of the center bin CpG with the TF bound or nucleosome bound fraction.
* For TFBS the JASPAR 2018 database is used and curated using ChIP-seq data.
* For more information on SMF analysis of TFBS, see [Kleinendorst and Barzaghi et al., 2022]() and [Soenmezer et al., 2022]().

# Covered analyses
Scripts for the following analyses are provided:
* Single molecule methylation call - CA strategy 
* blabla




# Preamble: preprocessing SMF data
To ensure compatibility with our downstream tools, we recommend aligning sequencing reads using the QuasR function [qAlign](https://www.rdocumentation.org/packages/QuasR/versions/1.12.0/topics/qAlign) as follows
```r
cl = makeCluster(40)
prj = QuasR::qAlign(sampleFile = sampleFile,
              genome = genome,
              aligner = "Rbowtie",
              projectName = "prj", 
              paired = "fr",
              bisulfite = "undir", 
              alignmentParameter = "-e 70 -X 1000 -k 2 --best -strata",
              alignmentsDir = "./", 
              cacheDir = tempdir(),
              clObj = cl)
```
For more details on how to structure the **sampleFile** argument we refer to the [qAlign](https://www.rdocumentation.org/packages/QuasR/versions/1.12.0/topics/qAlign) documentation.
For more details on SMF data preprocessing we refer to the computational steps of our SMF methods manuscript [Kleinendorst and Barzaghi et al., 2022]().

# Installation
???????

# Descriptions of the analyses
More detailed descriptions of the individual analyses can be found in the STAR methods section of the publication [Kreibich et al., 2022, bioRxiv](https://doi.org/10.1101/2022.05.19.492653).
