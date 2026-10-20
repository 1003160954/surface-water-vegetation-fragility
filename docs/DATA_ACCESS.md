# Data access and provenance

Large raw exports are retained by the authors outside this repository. See reviewer_evidence/raw_input_manifest.json for the analysis export's file name, bytes, exact SHA256, row counts and filter. It is a processed panel export rather than raw satellite observations. Original EVI/LAI/thermal/ERA5 collection details are not all recovered, so product families must not be substituted for verified IDs.

Derived national corrected flags, harmonized covariates, event numerators, model outputs and figure-source tables are included. They permit arithmetic verification and result visualization. Fitting models and reconstructing every map still needs the original exported panel. No public access to restricted Drive files is granted by this repository, and no credentials or large source parquet is included.

The source-only notebooks are recovered code archives without historical outputs. Current reproducibility includes exact selected downstream results; it does not mean the entire original upstream extraction can be rebuilt independently. Dataset catalog links and attribution requirements are recorded in docs/data_sources.csv and README.md. WorldClim background covariates are different from an independent event-anomaly reference.
