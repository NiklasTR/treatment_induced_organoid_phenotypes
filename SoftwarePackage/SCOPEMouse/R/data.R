#' Viability Classifier Accuracies
#'
#' The full accuracy matrix of classifiers for individual organoid lines (rows)
#' applied to features of all other lines (columns)
#'
#' @format A matrix with 4 rows and columns
#' \describe{
#'   \item{rows}{Viability classifiers trained on each organoid line}
#'   \item{columns}{Validation data for each organoid line}
#' }
"acc_matrix"

#' Viability Classifier Accuracies for Two-Channel Classifiers
#'
#' A list of the full accuracy matrices of classifiers for individual organoid
#' lines (rows) applied to features of all other lines (columns). The list has
#' four entries, one for each reduced classifier: only actin, only cell nuclei,
#' only dead cells, or actin + cell nuclei.
#'
#' @format A list with 4 entries. Each entry is a 4 x 4 matrix.
"acc_matrix_reduced"

#' Drug Target Annotations
#'
#' Pathways and targets for all drugs in the screen.
#'
#' @format Data frame with 1699 rows and 3 columns
#' \describe{
#'   \item{Drug}{}
#'   \item{Pathways}{}
#'   \item{Targets}{}
#' }
"drug_annotations"

#' Viabilities of wells
#'
#' Classification results for wells.
#'
#' @format Data frame with 18400 rows and 12 columns
#' \describe{
#'   \item{Product.Name}{}
#'   \item{Concentration}{}
#'   \item{Plate.ID}{}
#'   \item{Well.ID}{}
#'   \item{Num.Objects}{}
#'   \item{Percent.Dead}{Percent of organoids that are dead in the well}
#'   \item{Percent.Live}{Percent of organoids that are alive in the well}
#'   \item{Mean.Certainty.Live}{The mean classification accuracy of all organoids classified as alive}
#'   \item{Mean.Certainty.Dead}{The mean classification accuracy of all organoids classified as dead}
#'   \item{Line}{}
#'   \item{Replicate}{}
#'   \item{Layout}{}
#' }
"mortality"

#' Viabilities of wells
#'
#' A list of classification results for wells for all reduced-channel
#' classifiers. The list has four entries, one for each reduced classifier:
#' only actin, only cell nuclei, only dead cells, or actin + cell nuclei.
#'
#' @format A list with 4 entries. Each entry is a data frame with 18400 rows
#' and 12 columns.
#' \describe{
#'   \item{Product.Name}{}
#'   \item{Concentration}{}
#'   \item{Plate.ID}{}
#'   \item{Well.ID}{}
#'   \item{Num.Objects}{}
#'   \item{Percent.Dead}{Percent of organoids that are dead in the well}
#'   \item{Percent.Live}{Percent of organoids that are alive in the well}
#'   \item{Mean.Certainty.Live}{The mean classification accuracy of all organoids classified as alive}
#'   \item{Mean.Certainty.Dead}{The mean classification accuracy of all organoids classified as dead}
#'   \item{Line}{}
#'   \item{Replicate}{}
#'   \item{Layout}{}
#' }
"mortality_reduced"

#' ROC curves for viability classifiers
#'
#' ROC curves for viability classifiers for all combinations of classifiers
#' applied to data of other lines.
#'
#' @format Data frame with 5 columns
#' \describe{
#'   \item{Threshold}{The classification threshold being tested. Sometimes, a
#'   threshold of 2.0 is selected by python. This is only a technical measure
#'   to ensure that the ROC curve goes from (0, 0) to (1, 1).}
#'   \item{FalsePosRate}{False positive rate at given threshold}
#'   \item{TruePosRate}{True positive rate at given threshold}
#'   \item{clf_line}{The source line of the classifier being tested}
#'   \item{data_line}{The source line for the data being tested on}
#' }
"roc_data"

#' ROC curve AUCs
#'
#' AUC values for ROC curves found in 'roc_data'
#'
#' @format Data frame with 3 columns
#' \describe{
#'   \item{clf_line}{The source line of the classifier being tested}
#'   \item{data_line}{The source line for the data being tested on}
#'   \item{AUC}{Area under receiver operating curve}
#' }
"roc_aucs"

#' Drug-induced phenotype profiles
#'
#' Drug-induced phenotype vector as calculated by the SVMs. Each row in the
#' matrix corresponds to one treatment. A treatment consists of a line, drug,
#' and concentration (when applicable). The treatment ID can be read from the
#' row names with the format 'LINE.DRUG__CONCENTRATION' (for Clinical Cancer
#' Panel drugs) or 'LINE.DRUG' (for KiStem drugs). 'drug_effect_metadata' also
#' has this information in a more easily manageable table format.
#'
#' Profiles are averaged over 10 cross validation iterations.
#'
#' @format Matrix with 6796 rows and 25 columns
"drug_effect_profiles"

#' Drug-induced phenotype profile metadata
#'
#' Metadata belonging to 'drug_effect_profiles'. Row X in the metadata
#' corresponds to row X in the profiles (rownames are identical)
#'
#' @format Data frame with 6796 rows and 5 columns
#' \describe{
#'   \item{AUC_Mean}{The mean SVM AUROC, i.e. classification performance,
#'   averaged over 10 cross validation iterations.}
#'   \item{AUC_Std}{The standard deviation of the SVM AUROC over 10 cross
#'   validation iterations.}
#'   \item{Distance}{The euclidean distance between the median features of
#'   negative controls and treated drug (after PCA-transformation).}
#'   \item{Line}{}
#'   \item{Drug}{}
#' }
"drug_effect_metadata"

#' Drug-induced phenotype profiles when computed replicate-wise
#'
#' Drug-induced phenotype vector as calculated by the SVMs. In this version,
#' replicates are treated as separate treatments. Each row in the
#' matrix corresponds to one treatment. A treatment consists of a line, drug,
#' and concentration (when applicable). The treatment ID can be read from the
#' row names with the format 'LINE.DRUG__CONCENTRATION' (for Clinical Cancer
#' Panel drugs) or 'LINE.DRUG' (for KiStem drugs). 'drug_effect_metadata' also
#' has this information in a more easily manageable table format.
#'
#' Profiles are averaged over 10 cross validation iterations.
#'
#' @format Matrix with 13592 rows and 25 columns
"drug_effect_rep_profiles"

#' Drug-induced phenotype profile metadata
#'
#' Metadata belonging to 'drug_effect_rep_profiles'. Row X in the metadata
#' corresponds to row X in the profiles (rownames are identical)
#'
#' @format Data frame with 13592 rows and 5 columns
#' \describe{
#'   \item{AUC_Mean}{The mean SVM AUROC, i.e. classification performance,
#'   averaged over 10 cross validation iterations.}
#'   \item{AUC_Std}{The standard deviation of the SVM AUROC over 10 cross
#'   validation iterations.}
#'   \item{Distance}{The euclidean distance between the median features of
#'   negative controls and treated drug (after PCA-transformation).}
#'   \item{Line}{}
#'   \item{Drug}{}
#' }
"drug_effect_rep_metadata"

#' Well-averaged features
#'
#' Well-averages of processed single-organoid features. Note that this still
#' contains all features, even those that were too constant, i.e. resulted
#' in a value of NA/NaN after preprocessing. Row names indicate plate and
#' well ID.
#'
#' @format Matrix with 18416 rows and 3143 columns
"well_features"

#' Well-averaged features metadata
#'
#' Metadata belonging to 'well_features'.
#'
#' @format Matrix with 18416 rows and 6 columns
"well_metadata"
