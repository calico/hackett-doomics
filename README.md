# Hackett et al. 2023 - BioRxiv

This study was a collaboration with Calico Life Sciences and Gary Churchill aimed to identifying circulating molecules which are correlated with age and/or lifespan. To do this, plasma was profiled at three ages in each of 110 mice using proteomics, metabolomics, and lipidomics. 

This repository hosts the codebase used to carryout the studies major analyses and also serves as a landing page for other resources (apps and data) associated with the manuscript.

## Shiny App

Any of the ~2,300 molecular features described in this study can be visualized using our [Shiny App](http://public-rstudio-connect.calicolabs.com/doomics_results_browser/). This application allow you to look at how a molecule varies with age, fraction of life lived, lifespan, etc, and to stratify results by possible confounding effects such as experimental batches.

## Repository structure

Analyses are divided into two major phases:

- Bioinformatics: Code and notebooks which generated the data present in "Google Cloud Storage (GCS)" but which require infrastructure (e.g., a SLURM cluster for bootstrapping), or logistics (e.g., interfacing between data in googleDrive, GCS, Cromwell), which cannot be easily reproduced by end-users. For clarity, notebook have been rendered and can be seen below:
  - 1_wdl_processes: metabolomics and lipidomics informatics run through Cromwell
    - Metabolomics: set1/2 [**( + )**](http://public-rstudio-connect.calicolabs.com/doomics_met_pos_12/), [**( - )**](http://public-rstudio-connect.calicolabs.com/doomics_met_neg_12/), set3 [**( + )**](http://public-rstudio-connect.calicolabs.com/doomics_met_pos_3/), [**( - )**](http://public-rstudio-connect.calicolabs.com/doomics_met_neg_3/)
    - Lipidomics: set1/2 [**( + )**](http://public-rstudio-connect.calicolabs.com/doomics_lipid_pos_12/), [**( - )**](http://public-rstudio-connect.calicolabs.com/doomics_lipid_neg_12/), set3 [**( + )**](http://public-rstudio-connect.calicolabs.com/doomics_lipid_pos_3/), [**( - )**](http://public-rstudio-connect.calicolabs.com/doomics_lipid_neg_3/)
  - 2_non_wdl_processes: manually run notebooks which aggregate proteomics results and perform differential expression (with bootstrapping in SLURM):
    - [Proteomics](http://public-rstudio-connect.calicolabs.com/doomics_proteomics/)
    - [Differential Expression](http://public-rstudio-connect.calicolabs.com/doomics_diffex/)
  - 3_interactive: organize results and upload them to GCS to support major_analyses and the Shiny app
    - [GCS Setup](http://public-rstudio-connect.calicolabs.com/gcs_setup/)

- Major Analyses: Reproducible notebooks reproducing the vast majority of figures and tables starting from data which end-users can download from Google Cloud Platform (GCP). These notebooks roughly mirror the structure of the manuscript with individual notebooks indicating the section of the manuscript they support:

  1. [Figure1 - Study Design](http://public-rstudio-connect.calicolabs.com/doomics_figure1/)
  2. [Figure2 - Exploratory Data Analysis](http://public-rstudio-connect.calicolabs.com/doomics_figure2/)
      a. [Figure2 - Genetics](http://public-rstudio-connect.calicolabs.com/doomics_figure2_genetics/)
      b. [Figure2 - Batch Effects](http://public-rstudio-connect.calicolabs.com/doomics_figure2_batches/)
  3. [Figure3 - Statistics and Functional Enrichments](http://public-rstudio-connect.calicolabs.com/doomics_figure3/)
      a. [Figure3 - Power Analysis](http://public-rstudio-connect.calicolabs.com/doomics_figure3_power/)
  4. [Figure4 - Aging Archetypes](http://public-rstudio-connect.calicolabs.com/doomics_figure4/)
      a. [Figure3 - Volcano Plots](http://public-rstudio-connect.calicolabs.com/doomics_figure4_volcanos/)
  5. [Figure5 - Molecular Architecture of Longevity](http://public-rstudio-connect.calicolabs.com/doomics_figure5/)
  6. [Figure6 - Mechanisms](http://public-rstudio-connect.calicolabs.com/doomics_figure6/)

## How to run analyses

1. Set *config.json*:
  - cache_dir: set this to a location on your computer where you want raw data and cached intermediate files to be stored.
  - update_figure: this should be false (this option was just used to prepare the manuscript)
2. Load *hackett-doomics.Rproj* in Rstudio. (this will set your working directory to the project directory to normalize paths, and set the project R environment to use pre-set package versions via "renv").
3. Run *renv::restore()* to create install the packages required for all notebooks. This environment currently expected R **4.3.1**.
4. Running any of the notebooks in "major_analyses" will download and cache the studies intermediate data files.
5. major_analyses can be run in order to generate nearly all figures and tables in the manuscript. To do this run figure_1.qmd > figure_2.qmd > figure_3.qmd, etc. Analyses with a suffix, e.g., figure_2_batch.qmd are meant to be run after the corresponding primary file (in this case figure_2.qmd). Rendered versions of all notebooks can be seen above.

