# Congruent SSEs
 The reports and R scripts for demonstrating model congruence in state-dependent diversification models (SSEs).

## Preprint:

Tarasov, S. and Uyeda, J., 2024. Nonidentifiability of state-dependent diversification models (SSEs) is ubiquitous but not problematic for phylogenetics. *bioRxiv*, [doi: https://doi.org/10.1101/2022.07.04.498736](https://www.biorxiv.org/content/10.1101/2022.07.04.498736v3)



## Quick start guide


- `Report` is the main directory to see the reports for all simulation performed by the study. It contains the following sections:
  - 01-Validation.Rmd
  - 02-MiSSE.Rmd
  - 03-CID.Rmd
  - 04-EHE.Rmd
  - 05-Semi_Congruent.Rmd
  - 06-Empirical.Rmd 
  - 07-Sampling_Fraction.Rmd
  - 08-Sim_Data.Rmd
  - 09-Topo_bridges.Rmd
- `Rmd` contains the code to reprodice the reports and simulation.
- `R` contains necessary R functions.
  - `/R/hiclasse/src/` includes C code for HiClaSSE ODEs.
- `Mathematica` includes Mathematica Notebook that creates symbolic ODEs for HiClaSSE and can be used to generate C code.



 <p align="left">
  <img src="https://github.com/sergeitarasov/Congruent-SSE-CTMC/blob/main/Fig_class.png" width="600" title="hover text">
</p>  
