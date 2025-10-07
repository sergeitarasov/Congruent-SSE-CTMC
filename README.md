
# Congruent-SSE-CTMC

This repository contains R scripts and notebooks associated with the manuscript *Universal Congruence in State-Dependent Diversification Models (SSEs): Resolving Unidentifiability and Explaining High Error Rates*.  
It demonstrates model congruence in state-dependent diversification (SSE) models and provides the code and data necessary to reproduce the main results.  

A Zenodo-archived snapshot of this repository is available at:  
https://doi.org/10.5281/zenodo.17256718  

---

## 1. System Requirements

- **Operating System:** Tested on macOS (15.7), but should work equally on Linux and Windows.  
- **R version:** ≥ 4.5.1  
- **Dependencies:** No dependencies other than R and the associated R packages.  
- **Hardware:** No special hardware required. Typical analyses run on a standard laptop/desktop with 8 GB RAM.  

---

## 2. Installation Guide

1. Install R (≥ 4.5.1) from [CRAN](https://cran.r-project.org/).  
2. Clone or download this repository:  
   ```bash
   git clone https://github.com/sergeitarasov/Congruent-SSE-CTMC.git
   ```
   Or download the ZIP archive and extract it.  
3. Navigate to the `Rmd/` folder to run the R notebooks.  

**Typical install time:** ~1 minute on a standard desktop computer.   

---

## 3. Demo

A quick demo can be run using `Rmd/02-BiSSE.Rmd`.  

**To run the demo:**  
1. Open `Rmd/02-BiSSE.Rmd` in RStudio.  
2. Execute the code blocks.  

**Expected output:**  
- Model fitting of BiSSE, congruent and non-congruent models.  
- Plots of likelihood values.  

**Expected run time:** < 5 minutes on a standard desktop computer.  

---

## 4. Instructions for Use

- Execute all Rmd files in the `Rmd/` folder; they contain the necessary instructions and code to reproduce the analyses.    

---

## 5. Content

- `Report` is the main directory to see the html reports for all simulation performed by the study. It contains the following sections:
  - 01-Validation.html
  - 02-BiSSE.html
  - 03-CID.html
  - 04-EHE.html
  - 05-Semi_Congruent.html
  - 06-Empirical.html
  - 07-Sampling_Fraction.html
  - 08-Sim_Data.html
  - 09-Topo_bridges.html
- `Rmd` contains the code to reprodice the reports and simulation.
- `R` contains necessary R functions.
  - `/R/hiclasse/src/` includes C code for HiClaSSE ODEs.
- `Mathematica` includes Mathematica Notebook that creates symbolic ODEs for HiClaSSE and can be used to generate C code (it is not necessary for reproducing the results).



 <p align="left">
  <img src="https://github.com/sergeitarasov/Congruent-SSE-CTMC/blob/main/Fig_class.png" width="600" title="hover text">
</p>  
