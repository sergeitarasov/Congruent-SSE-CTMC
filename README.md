# Congruent SSEs
 The reports and R scripts for demonstrating model congruence in state-dependent diversification models (SSEs).

## Preprint:

Tarasov, S. and Uyeda, J., 2024. Ubiquitous nonidentifiability of state-dependent diversification models (SSEs) is not problematic for phylogenetics. *bioRxiv*, [doi: https://doi.org/10.1101/2022.07.04.498736](https://www.biorxiv.org/content/10.1101/2022.07.04.498736v3)



## Quick start guide


- `Report.zip` is the main file to see the reports for all simulation performed by the study. This file is the zipped version of the `Report` folder. You need to download it, unzip it and click on index.html to browse through the reports. It contains the following sections:
  - Setting-up HiClaSSE
  - Validation of HiClaSSE
  - MiSSE Congruence
  - Congruence between independent and dependent SSEs
  - Equal Rate Hidden Expansion (EHE)
  - Semi-congruent Behavior
  - Congruence under varying Sampling Fraction
- `Rmd` contains the code to reprodice the reports and simulation.
- `R` contains necessary R functions.
  - `/R/hiclasse/src/` includes C code for HiClaSSE ODEs.
- `Mathematica` includes Mathematica Notebook that creates symbolic ODEs for HiClaSSE and can be used to generate C code.



 <p align="left">
  <img src="https://github.com/sergeitarasov/Congruent-SSE-CTMC/blob/main/Fig_class.png" width="600" title="hover text">
</p>  
