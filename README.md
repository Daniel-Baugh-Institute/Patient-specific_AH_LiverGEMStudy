# Patient-specific alcohol-associated hepatitis liver genome-scale metabolic modeling study

This repository contains 10 files.
- **CalculateNetEnzymaticLevel.R**: R script for calculating the Net Enzymatic Level (NEL) for each reaction in the Human1 generic GEM
- **Hsapiens_cellMedium.csv**: CSV file containing the Homo Sapien cell medium information necessary for running the Pheflux algorithm
- **LiverMetabolomics.SevereAH.vs.HealthyControls.xlsx**: Excel file containing the metabolites that are significantly different between Severe AH and Healthy individuals
- **PhefluxInputData_Example.csv**: CSV file showing an example input file neccessary for running the Pheflux algorithm
- **Analysis_and_FigureReplication.R**: R script for replicating the figures and performing the analysis in the manuscript
- **Run_PheFlux.ipynb**: Jupyter notebook file that contains the code for running Pheflux
- **SaveFluxFiles_PhefluxOutput.R**: R script for saving the Pheflux output files as one large matrix
- **SaveGeneExpFiles_forPheflux.R**: R script for saving the gene expression files necessary for running the Pheflux algorithm
- **RstudioInformation.txt**: TXT file with information on the RStudio working environment used for analysis and figure replication in the manuscript
- **PythonInformation.txt**: TXT file with information on the Python working environment used for running the Pheflux algorithm

- # Dependencies
- [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) model (ver. 1.6.0)
- [PheFlux](https://github.com/mrivas/pheflux) algorithm (ver. 1.0.0)

**ref:** Manchel et al. Genome-scale patient-specific modeling identifies combinatorial intervention via *Hkdc1* and *Pkm* to reverse glucose use metabolic reprogramming in alcohol-associated hepatitis. In Preparation.
