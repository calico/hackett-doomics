## Execution

Analyses are organized and run using Cromwell. Download the most recent version, e.g.,:

    wget https://github.com/broadinstitute/cromwell/releases/download/52/cromwell-52.jar

Setup your environment:
- Add to .bash_profile and source:

     doomics_repo_path=<<path_to_repo>>

And then execute all analyses:

    java -jar ~/apps/cromwell-52.jar run ${doomics_repo_path}/bioinformatics/1_wdl_processes/plasma_omics.wdl
            
            
## Results

A versioned informatics pipeline, analyses and supporting data files is stored as a single Cromwell directory. Each subdirectory corresponds to a .Rmd notebook that was run to create a .html, .png versions of the figures from the notebook (out_figures) and any outputs of the analyses (out_files). 

This includes:

- **small-molecule-xxx** an standardized analysis of a single metabolomics or lipidomics dataset
- **differential_abundances** - testing whether molecular features change based on age, lifespan, ...
    - out_files
        - model_data.Rds - all measurements of all data modalities in tidy format including sample metadata.
        - lm_fits.Rds - fitted values (and residuals) of each model. This is useful for EDA.
        - target_model_parameters.Rds/tsv - parameter estimates from regression models
        - model_significance.Rds - target_model_parameters with FDR and feature metadata