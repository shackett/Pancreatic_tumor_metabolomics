# Pancreatic_tumor_metabolomics
Scripts used to compare metabolite concentrations between paired tumor-benign human pancreatic cancer samples. The main results of this analysis are presented in Figure 1 in Kamphorst et al. 2015.

Jurre Kamphorst, Michel Nofal, Cosimo Commisso, Sean R. Hackett, Wenyun Lu, Elda Grabocka, George Miller, Jeffrey Drebin, Matthew Vander Heiden, Dafna Bar-Sagi, Craig Thompson, Josh Rabinowitz. [Human pancreatic cancer tumors are nutrient poor and the tumor cells actively scavenge extracellular protein](http://cancerres.aacrjournals.org/content/75/3/544.long). Cancer Research, 75, 2015.

------

cancer_analysis_stream.R : Analysis pipeline
- Load data and separate information into metabolite-specific information, sample-specific information and a metabolite abundance matrix
- Exclude patients whose tumors are not PDAC (Pancreatic ductal adenocarcinoma)
- Using log-ion counts, each sample (a single mass-spec run on a single instrument) is corrected for loading variability using a robust median polish
- For each patient x tumor/benign, a point estimate of log-ion counts is found by taking the mean over replicates
- For exploratory analysis every peak, on every instrument is retained
- For non-exploratory analysis, if a 1-1 match based on m/z and RT between a metabolite and a peak could not be found, this metabolite was excluded from analysis (on that instrument).  For metabolites that are measured on multiple instruments, the signals from each instrument are combined to generate a concensus
- For each metabolite and patient, the fold-change between their tumor metabolite concentration and benign metabolite concentration is found.
- For each metabolite, tumor-benign differences are bootstrapped to determine whether the average fold-change of a metabolite's concentration differs significantly from the expectation of zero.
- P-values were calculating assuming a two-tailed test and these p-values were FDR corrected to generate q-values
- Discoveries at an FDR of 0.05 were chosen for further analysis


cancer_lib.R : Helper functions
- Input reformating: the input data frame is converted into 3 pieces:
  - TNmeta: Data matrix (n x m): raw ion counts
  - meta.data: Metabolite information data.frame (n rows): metabolite name, instrument, etc.
  - header.data: Sample information data.frame (m rows): individual, date, tumor/benign
- QC functions:
  - test the distribution of T-N residuals
  - determine whether *ex vivo* time affects metabolite abundance to test metabolite stability
- Plotting: boxplots and heatmaps
