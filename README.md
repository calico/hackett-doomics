# Hackett et al. 2023 - BioRxiv

This study was a collaboration with Calico Life Sciences and Gary Churchill aimed to identifying circulating molecules which are correlated with age and/or lifespan. To do this, plasma was profiled at three agees in each of 110 mice using proteomics, metabolomics, and lipidomics. 

This repository hosts the codebase used to carryout the studies major analyses and also serves as a landing page for other resources (apps and data) associated with the manuscript.

## Repository structure

Analyses are divided into two major phases:

- Major Analyses: Reproducible notebooks reproducing the vast majority of figures and tables starting from data which end-users can download from Google Cloud Platform (GCP). 
- Bioinformatics: Code and notebooks which generated the data present in "Google Cloud Storage (GCS)" but which require infrastructure (e.g., a SLURM cluster for bootstrapping), or logistics (e.g., interfacing between data in googleDrive, GCS, Cromwell), which cannot be easily reproduced by end-users.

## How to run analyses

1. Set config.json:
  - cache_dir: set this to a location on your computer where you want raw data and cached intermediate files to be stored
2. Load hackett-doomics.Rproj in Rstudio. (this will set your working directory to the project directory to normalize paths, and set the project R environment to use pre-set package versions via "renv").
3. Running any of the notebooks in "major_analyses" will download and cache the studies intermediate data files.
4. major_analyses can be run in order to generate nearly all figures and tables in the manuscript. To do this run figure_1.qmd > figure_2.qmd > figure_3.qmd, etc. Analyses with a suffix, e.g., figure_2_batch.qmd are meant to be run after the corresponding primary file (in this case figure_2.qmd). Rendered versions of all notebooks can be seen above.

## Dev Zone [meant for the studies authors]




