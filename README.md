# Patient-specific alcohol-associated hepatitis liver genome-scale metabolic modeling study

This repository contains 11 files and 1 directory. <br /><br />
**Files:**
- **CalculateNetEnzymaticLevel.R**: R script for calculating the Net Enzymatic Level (NEL) for each reaction in the Human1 generic GEM
- **Hsapiens_cellMedium.csv**: CSV file containing the Homo Sapien cell medium information necessary for running the Pheflux algorithm
- **Human-GEM.xml**: XML file of the Human1 generic metabolic network
- **LiverMetabolomics.SevereAH.vs.HealthyControls.xlsx**: Excel file containing the metabolites that are significantly different between Severe AH and Healthy individuals
- **PhefluxInputData_Example.csv**: CSV file showing an example input file neccessary for running the Pheflux algorithm
- **Analysis_and_FigureReplication.R**: R script for replicating the figures and performing the analysis in the manuscript
- **Run_PheFlux.ipynb**: Jupyter notebook file that contains the code for running Pheflux
- **SaveFluxFiles_PhefluxOutput.R**: R script for saving the Pheflux output files as one large matrix
- **SaveGeneExpFiles_forPheflux.R**: R script for saving the gene expression files necessary for running the Pheflux algorithm
- **RstudioInformation.txt**: TXT file with information on the RStudio working environment used for analysis and figure replication in the **Analysis_and_FigureReplication.R** script
- **PythonInformation.txt**: TXT file with information on the Python working environment used for running the Pheflux algorithm in the **Run_PheFlux.ipynb** script <br />

**Directory:**
- **FilesforFigureReplication**: Directory containing all necessary files for analysis and figure replication using the **Analysis_and_FigureReplication.R** script

# Dependencies
- RStudio (ver. 4.2.2)
- Python (ver. 3.9.7)
- [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) model (ver. 1.6.0)
- [PheFlux](https://github.com/mrivas/pheflux) algorithm (ver. 1.0.0)

# Replication:
1. Install Python (ver. 3.9.7) and all Python dependencies in "PythonInformation.txt".
2. Install RStudio (ver. 4.2.2) and all RStudio dependencies in "RstudioInformation.txt".

**Patient-specific metabolic flux analysis:** 

1. Run "SaveGeneExpFiles_forPheflux.R" script to write patient-specific gene expression data to csv files for generating patient-specific metabolic fluxes using Pheflux.
2. Run "Run_PheFlux.ipynb" script to generate patient-specific metabolic fluxes using Pheflux.
3. Run "SaveFluxFiles_PhefluxOutput.R" to save all patient-specific metabolic fluxes to one large matrix.

**Downstream Analysis and Figure Replication:**

Run "CalculateNetEnzymaticLevel.R" and "Analysis_and_FigureReplication.R" using files in the "FilesforFigureReplication" directory to replicate downstream analysis results and figures from the manuscript.

#
**ref:** Manchel et al., Genome-scale patient-specific modeling identifies combinatorial intervention via *Hkdc1* and *Pkm* to reverse glucose use metabolic reprogramming in alcohol-associated hepatitis. In Review. Communications Biology.
