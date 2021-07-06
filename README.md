promise
==============================

Human Cancer Organoid Profiling Project at the German Cancer Research Center


Hackathon agenda: 
* manuscript
    * last version
    * comments
    * upcoming version
* overview github repo
    * separate existence of SCOPE package with rest of repository (especially relevant for package data within SCOPE)
    * ground truth file for organoid classification
    * RNA expression data harmonization


Issues I needed to fix
* document the analysis workflow
* run code with latest freetype version to perform ggrastr plots -> I add a subset of the data for local processing

Note for contributors
------------
the most up-to-date branch is "niklas". Please merge your contributions to this branch.
I plan to merge niklas and master in 08/21

To-Do
------------
File transfer:
* ctg data


Gene expression
* data were generated via the DKFZ GPCF
* raw data is stored in .CEL format under data/raw/expression, data relating to the latest manuscript is separated under microarray_manuscript
* a normalized file, also deposited at GEO, is generated
    
    
    
    
Analysis Pipeline
------------


raw.hdf -- FeatureAnalysis -> pca.hdf & raw.hdf
pca.hdf -- UMAP ->


Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


## data/raw/PROMISE (imaging)
* hdf5dir - contains raw well-level hdf5 imaging data of each well (e.g. per object: 4 regions, 3 channels, 16 z-levels), ca. 1GB/well
* hdf5projection - contains hdf5 objects with image projection information on a well level, ca. 50MB/well
* hdf5validation - text files containing error reports in case checksum are not matching between hdf5 files after file transfer (data was originally moved after local projection to a remote compute environment and storage)
* htmldir - visualization of plate projections
* segmentation - image segmentation segmentation masks on well-level
* segmentationhtmldir - visualization of segmentation masks
* features - storage of feature extraction results on a plate level with 3 elements
	* wells/ - directory with well-level feature data
	* features - a plate level aggregation of well-level features
	* processedFeatures - these data have been pre-processed, including the following steps
		* removing objects at the well boundary
		* removing objects that are out of focus
		* removing objects with unexpected size

## data/interim (created by ML tools/ scripts)
the directory contains 3 larger ML projects that were started for particular purposes within the manuscript
* **line_differences** - a collection of common features across all organoid lines, followed by PCA for Unsupervised Learning and EDA
	* Features_human_all_drugs.h5 - 27GB large file containing all features across lines
	* results/ReducedFeatures_all_drugs_human.h5 - 27GB large file containing all shared features across lines
	* results/ReducedFeaturesPCA_all_drugs_human.h5 - 5G large file containing results of a incremental PCA with 25 preserved components across all analyzed lines
* **drug_effects** - a project to estimate the effect of each drug treatment on every individual organoid line. Common features across lines are identified, features for each line are scaled to their DMSO control and PCA is being performed.
	* TransformedFeatures_human_25components.h5 - 5G large file containing results of a incremental PCA with 25 preserved components across all analyzed lines, the **important difference to ReducedFeaturesPCA_all_drugs_human.h5** is that input features were scaled to each line's respective DMSO control.
	* incrementalPCA - pickle of PCA model
	* lines/ - containing the linear models for each line and drug with their respective AUROC and pvalue
* **organoid_viability**
	* classifiers/ - model checkpoint of random forest
	* diagnostics/ 
	* results/ - containing csv files with viability estimates

## data/processed (data for notebook processing)
* umap_ ..
	* umap_absolute_all_drugs_tidy.Rds - all included organoid objects in UMAP projection with metadata
	* umap_absolute_all_drugs_sampled.Rds - a subset of the UMAP object above, that represents 5% of the original corpus, ca. 300k objects


<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
