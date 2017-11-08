---
title: "Screen Quality Control (Mouse)"
author: "Jan Sauer"
date: "22.10.2017"
output: 
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
---

This vignette performs quality control of the mouse screen.



# Data preprocessing

The data preprocessing consists of the following steps:

* Separate objects into "organoids" (area >= 2500 pixels) and "shrapnel" (area <= 2500 pixels).
* Calculate the summary features for each of these categories of objects within a well (median / mad).
  * I decided for this summary statistic as some features may only have features within a certain range, e.g. [0, 1] for the asphericity. Even under the assumption of normally distributed variation within a well, this limitation would cause the mean to skew the summary statistic towards the middle of the interval.
* Subtract from each feature the median of the DMSO wells and divide by their standard deviation on each plate to eliminate batch effects.
* Combine plates belonging to a single cell line and calculate the z-score for each feature.

This is done in an external python script. The data and all intermediate steps are stored in the corresponding hdf5 files in the PROMISE features directory.


```r
cell_lines = unique(substr(list.files(feature_dir, pattern = "M001"), 1, 7))
features = setNames(
  object = vector(mode="list", length=length(cell_lines)), 
  nm = cell_lines)
features_metadata = setNames(
  object = vector(mode="list", length=length(cell_lines)), 
  nm = cell_lines)
for(cell_line in cell_lines) {
  cl_feature_fn = file.path(feature_dir, sprintf("%s_averaged_features_%s.h5", cell_line, feature_type))
  features[[cell_line]] = data.frame(h5read(cl_feature_fn, hdf5_dataset))
  colnames(features[[cell_line]]) = h5read(cl_feature_fn, "feature_names")
  rownames(features[[cell_line]]) = h5read(cl_feature_fn, "well_names")
  
  features_metadata[[cell_line]] = data.frame(
    "REPLICATE" = h5read(cl_feature_fn, "replicates"),
    "CONCENTRATION" = h5read(cl_feature_fn, "concentrations"),
    "DRUG" = h5read(cl_feature_fn, "drugs"), 
    row.names = rownames(features[[cell_line]]))
}
features = do.call(rbind, features)
rownames(features) = sapply(strsplit(rownames(features), ".", fixed = TRUE), "[[", 2)
features_metadata = do.call(rbind, features_metadata)
rownames(features_metadata) = sapply(strsplit(rownames(features_metadata), ".", fixed = TRUE), "[[", 2)
```

For any number of reasons, some wells may not have been processed properly. In particular, if there is no information to the total number of organoids in a well then it can be safely discarded as this feature should always be computable.

```r
na_rows = rownames(features)[!is.finite(features$organoids_num.of.objects)]
na_rows_repid = paste0(substr(na_rows, 1, 7), "_", substr(na_rows, 12, 19))
all_wells_repid = paste0(substr(rownames(features), 1, 7), "_", substr(rownames(features), 12, 19))
features = features[!all_wells_repid %in% na_rows_repid, ]
features_metadata = features_metadata[!all_wells_repid %in% na_rows_repid, ]
```

Due to the normalization method, features may be marked as NaN or Infinity. This can happen in particular if a feature was constant across all DMSO wells of a plate. Since the wells are already averages over all organoids within the wells, I assume that a constant feature value across all DMSO wells indicates a bad feature. As I've already removed any bad wells, I now remove all features that do not have finite values in all remaining wells.

```r
nanfeatures = colSums(!is.finite(as.matrix(features)))
features = features[,nanfeatures == 0]
```

The features should now all be finite

```r
cat("Number of NaN features:", sum(!is.finite(as.matrix(features))))
```

```
## Number of NaN features: 0
```

Lastly, some of the wells were found to be entirely blurry and are effectively unusable. These have been previously classified with a random forest classifier. The blurry wells and their replicates are filtered out.

```r
blurry_wells = read.table("blurry_wells_mouse.txt", stringsAsFactors = FALSE)
blurry_wells_repid = paste0(substr(blurry_wells$V1, 1, 7), "_", substr(blurry_wells$V1, 12, 19))
all_wells_repid = paste0(substr(rownames(features), 1, 7), "_", substr(rownames(features), 12, 19))
wells_to_keep = !all_wells_repid %in% blurry_wells_repid
features = features[wells_to_keep,]
features_metadata = features_metadata[wells_to_keep,]
```

# Quality Control
## Cell Line Quality
Normally, one would look at the z-factor, a degree of separation between positive and negative controls, for each feature to determine the assay quality. As there are no known positive controls for the mouse screen, I instead look at the z-factor for all treatments compared to the negative controls (DMSO) for a few characteristic features, the number of organoids, the total biomass, and the average organoid size. I expect that at least some of the treatments should show considerable differences to the DMSO controls in these features.

The z-factor is defined as:

$$ Z^\prime = 1 - \frac{3 \cdot (\sigma_{treatment} + \sigma_{DMSO})}{|\mu_{treatment} - \mu_{DMSO}|}$$


```r
treatments = unique(features_metadata$DRUG)

all_biomass_zprime = list()
all_orgsize_zprime = list()
all_orgnum_zprime = list()

for(cell_line in cell_lines) {
  biomass_zprime = setNames(
    object = vector(mode="numeric", length=length(treatments)), 
    nm = treatments)
  orgsize_zprime = setNames(
    object = vector(mode="numeric", length=length(treatments)), 
    nm = treatments)
  orgnum_zprime = setNames(
    object = vector(mode="numeric", length=length(treatments)), 
    nm = treatments)
  dmso_treatment_biomass = features[
      features_metadata$DRUG == "DMSO" & 
        substr(rownames(features), 1, 7) == cell_line,
      "Total.Biomass"]
  dmso_treatment_orgsize = features[
      features_metadata$DRUG == "DMSO" & 
        substr(rownames(features), 1, 7) == cell_line,
      "organoids_x.0.s.area_expected"]
  dmso_treatment_orgnum = features[
      features_metadata$DRUG == "DMSO" & 
        substr(rownames(features), 1, 7) == cell_line,
      "organoids_num.of.objects"]
  for(treatment in treatments) {
    biomass_treatment_feat = features[
      features_metadata$DRUG == treatment & 
        substr(rownames(features), 1, 7) == cell_line, 
      "Total.Biomass"]
    orgsize_treatment_feat = features[
      features_metadata$DRUG == treatment & 
        substr(rownames(features), 1, 7) == cell_line, 
      "organoids_x.0.s.area_expected"]
    orgnum_treatment_feat = features[
      features_metadata$DRUG == treatment & 
        substr(rownames(features), 1, 7) == cell_line, 
      "organoids_num.of.objects"]
    
    biomass_zprime[[treatment]] = imageHTS::zprime(
      dmso_treatment_biomass,  biomass_treatment_feat, method = "robust")
    orgsize_zprime[[treatment]] = imageHTS::zprime(
      dmso_treatment_orgsize,  orgsize_treatment_feat, method = "robust")
    orgnum_zprime[[treatment]] = imageHTS::zprime(
      dmso_treatment_orgnum,  orgnum_treatment_feat, method = "robust")
  }
  
  all_biomass_zprime[[cell_line]] = biomass_zprime
  all_orgsize_zprime[[cell_line]] = orgsize_zprime
  all_orgnum_zprime[[cell_line]] = orgnum_zprime
}
```

```
## Error in loadNamespace(name): there is no package called 'imageHTS'
```

```r
all_biomass_zprime_df = do.call(cbind, all_biomass_zprime)
all_biomass_zprime_df = all_biomass_zprime_df[
  rowSums(!is.finite(as.matrix(all_biomass_zprime_df))) ==  0,]
```

```
## Error in array(x, c(length(x), 1L), if (!is.null(names(x))) list(names(x), : 'data' must be of a vector type, was 'NULL'
```

```r
all_biomass_zprime_df = apply(all_biomass_zprime_df, 2, sort)
```

```
## Error in apply(all_biomass_zprime_df, 2, sort): dim(X) must have a positive length
```

```r
all_biomass_zprime_df = melt(all_biomass_zprime_df)
```

```
## Error in names(object) <- nm: 'names' attribute [1] must be the same length as the vector [0]
```

```r
all_orgsize_zprime_df = do.call(cbind, all_orgsize_zprime)
all_orgsize_zprime_df = all_orgsize_zprime_df[
  rowSums(!is.finite(as.matrix(all_orgsize_zprime_df))) ==  0,]
```

```
## Error in array(x, c(length(x), 1L), if (!is.null(names(x))) list(names(x), : 'data' must be of a vector type, was 'NULL'
```

```r
all_orgsize_zprime_df = apply(all_orgsize_zprime_df, 2, sort)
```

```
## Error in apply(all_orgsize_zprime_df, 2, sort): dim(X) must have a positive length
```

```r
all_orgsize_zprime_df = melt(all_orgsize_zprime_df)
```

```
## Error in names(object) <- nm: 'names' attribute [1] must be the same length as the vector [0]
```

```r
all_orgnum_zprime_df = do.call(cbind, all_orgnum_zprime)
all_orgnum_zprime_df = all_orgnum_zprime_df[
  rowSums(!is.finite(as.matrix(all_orgnum_zprime_df))) ==  0,]
```

```
## Error in array(x, c(length(x), 1L), if (!is.null(names(x))) list(names(x), : 'data' must be of a vector type, was 'NULL'
```

```r
all_orgnum_zprime_df = apply(all_orgnum_zprime_df, 2, sort)
```

```
## Error in apply(all_orgnum_zprime_df, 2, sort): dim(X) must have a positive length
```

```r
all_orgnum_zprime_df = melt(all_orgnum_zprime_df)
```

```
## Error in names(object) <- nm: 'names' attribute [1] must be the same length as the vector [0]
```

```r
ggplotdf = data.frame(
  "Cell.Line" = all_biomass_zprime_df$Var2,
  "Total.Biomass" = all_biomass_zprime_df$value,
  "Num.Organoids" = all_orgnum_zprime_df$value,
  "Avg.Organoid.Size" = all_orgsize_zprime_df$value,
  "Treatment" = NA_real_)
```

```
## Error in data.frame(Cell.Line = all_biomass_zprime_df$Var2, Total.Biomass = all_biomass_zprime_df$value, : arguments imply differing number of rows: 0, 1
```

```r
ggplotdf = ggplotdf[
  is.finite(ggplotdf$Total.Biomass) & 
  is.finite(ggplotdf$Num.Organoids) &
  is.finite(ggplotdf$Avg.Organoid.Size), ]
```

```
## Error in eval(expr, envir, enclos): object 'ggplotdf' not found
```

```r
for(cell_line in cell_lines) {
  ggplotdf[ggplotdf$Cell.Line == cell_line, "Treatment"] = 
    seq_len(table(ggplotdf$Cell.Line)[cell_line])
}
```

```
## Error in table(ggplotdf$Cell.Line): object 'ggplotdf' not found
```


```r
ggplot(data = ggplotdf) +
  geom_line(aes(x = Treatment, y = Total.Biomass, color = Cell.Line)) + 
  coord_cartesian(ylim = c(-1, 1), xlim = c(1000, 1250)) + 
  ggtitle(label = "Z-Factor of the Total Biomass for the Treatments of each Cell Line")
```

```
## Error in ggplot(data = ggplotdf): object 'ggplotdf' not found
```


```r
ggplot(data = ggplotdf) +
  geom_line(aes(x = Treatment, y = Num.Organoids, color = Cell.Line)) + 
  coord_cartesian(ylim = c(-1, 1), xlim = c(1000, 1250)) + 
  ggtitle(label = "Z-Factor of the Number of Organoids for the Treatments of each Cell Line")
```

```
## Error in ggplot(data = ggplotdf): object 'ggplotdf' not found
```


```r
ggplot(data = ggplotdf) +
  geom_line(aes(x = Treatment, y = Avg.Organoid.Size, color = Cell.Line)) + 
  coord_cartesian(ylim = c(-1, 1), xlim = c(1000, 1250)) + 
  ggtitle(label = "Z-Factor of the Average Organoid Size for the Treatments of each Cell Line")
```

```
## Error in ggplot(data = ggplotdf): object 'ggplotdf' not found
```

While the quality of the screen is far from ideal, only cell line M001W01 gives reason for concern.

## Pruning Features

Features only make sense if they are sufficiently continuous. I want to discard features that are constant across the negative (DMSO) controls of a plate. A look at the 'minimum ratio' of unique values for each feature across plates shows that features are either continuous or constant across a plate. There is a reasonably steep transition between the two states so that the actual threshold chosen is not so important. I keep only features with more than 25% unique values across all plates.

```r
neg_ctrl = features[features_metadata$DRUG == "DMSO", ]
# This rounding ensures that there are no floating point arithmetic errors that are mistaken as unique values
neg_ctrl = round(neg_ctrl, 8)
neg_ctrl_cl = substr(rownames(neg_ctrl), 1, 14)
unique_vals = aggregate(neg_ctrl, list(neg_ctrl_cl), function(x) length(unique(x)) / length(x))
rownames(unique_vals) = unique_vals$Group.1
unique_vals$Group.1 = NULL
unique_vals_min = apply(unique_vals, 2, min)
unique_vals_min_df = data.frame(
  "Ratio.Unique.Vals" = sort(unique_vals_min, decreasing = TRUE), 
  "Features" = seq_along(unique_vals_min))
ggplot(data = unique_vals_min_df) + geom_point(aes(x = Features, y = Ratio.Unique.Vals)) + 
  ggtitle("Minimum Ratio of Unique Values for each Feature", 
          subtitle = "e.g. a 'minimum ratio' of 0.25 means all plates had a ratio of >= 25% of unique values for that feature") + 
  geom_hline(aes(yintercept = 0.25), color = "blue")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

```r
features_to_keep = names(which(unique_vals_min >= 0.25))
if(!"Total.Biomass" %in% features_to_keep)
  features_to_keep = c(features_to_keep, "Total.Biomass")
if(!"organoids_x.0.s.area_expected" %in% features_to_keep)
  features_to_keep = c(features_to_keep, "organoids_x.0.s.area_expected")
if(!"organoids_num.of.objects" %in% features_to_keep)
  features_to_keep = c(features_to_keep, "organoids_num.of.objects")
features = features[,features_to_keep]
```

## Correlation between replicates
### Cell line specfific
I look at all features now to see how well they correlate between replicates of the screen for each cell line individually


```r
cl_correlations = setNames(
  object = vector(mode = "list", length = length(cell_lines)), 
  nm = cell_lines
)
for(cell_line in cell_lines) {
  cl_features = features[substr(rownames(features), 1, 7) == cell_line, ]
  cl_features_metadata = features_metadata[substr(rownames(features_metadata), 1, 7) == cell_line, ]
  rep1 = cl_features[cl_features_metadata$REPLICATE == 1, ]
  rep2 = cl_features[cl_features_metadata$REPLICATE == 2, ]
  if(!identical(paste0(substr(rownames(rep1), 1, 7), "_", substr(rownames(rep1), 12, 19)), 
                paste0(substr(rownames(rep2), 1, 7), "_", substr(rownames(rep2), 12, 19)))) {
    warning(sprintf("Replicates for '%s' are not in the same order!", cell_line))
  }
  rep_correlations = setNames(
    object = vector(mode = "numeric", length = ncol(cl_features)), 
    nm = colnames(cl_features))
  for(feature in colnames(cl_features)) {
    rep_correlations[[feature]] = cor(
      rep1[[feature]], rep2[[feature]], 
      use = "pairwise.complete.obs")
  }
  rep_correlations = sort(rep_correlations, decreasing = TRUE)
  cl_correlations[[cell_line]] = data.frame(
    "Features" = seq_along(rep_correlations),
    "Feature.Name" = names(rep_correlations),
    "Correlations" = rep_correlations, 
    "Cell.Line" = cell_line
  )
}

cl_correlations = do.call(rbind, cl_correlations)
ggplot(data = cl_correlations) + geom_line(aes(x = Features, y = Correlations, color = Cell.Line)) + 
  ggtitle("Correlations between Replicates for each Individual Cell Line")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

Cell lines M001W01 and M001A03 both show inadequate correlations. This will become a problem during the feature selection, which requires that features correlate significantly between replicates.

### Screen-wide

```r
rep1 = features[features_metadata$REPLICATE == 1, ]
rep2 = features[features_metadata$REPLICATE == 2, ]
if(!identical(paste0(substr(rownames(rep1), 1, 7), "_", substr(rownames(rep1), 12, 19)), 
              paste0(substr(rownames(rep2), 1, 7), "_", substr(rownames(rep2), 12, 19)))) {
  warning("Replicates are not in the same order!")
}

screen_correlations = setNames(
  object = vector(mode = "numeric", length = ncol(features)), 
  nm = colnames(features))
for(feature in colnames(features)) {
  screen_correlations[[feature]] = cor(
    rep1[[feature]], rep2[[feature]], 
    use = "pairwise.complete.obs")
}

ggplotdf = data.frame(
  "Features" = seq_along(sort(screen_correlations)),
  "Correlation" = sort(screen_correlations, decreasing = TRUE))
ggplot(data = ggplotdf) + geom_line(aes(x = Features, y = Correlation)) + 
  ggtitle("Correlations between Replicates across all Cell Lines")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

Keep features that correlate between replicates within cell lines M001K02 and M001B03 (the other two exhibit insufficient correlation between replicates) with $\rho >= 0.6$. I also ensure that the number of organoids and the total biomass are guaranteed to be in the dataset (regardless of their correlation).

```r
cor_thresh = 0.6
cell_line_features = cl_correlations[
  cl_correlations$Correlations >= cor_thresh & 
  cl_correlations$Cell.Line %in% c("M001K02", "M001B03"), 
  "Feature.Name"]
cell_line_features = names(which(table(cell_line_features) == max(table(cell_line_features))))
features_to_keep = cell_line_features
# screen_features = names(screen_correlations[screen_correlations >= cor_thresh])
# features_to_keep = intersect(screen_features, cell_line_features)

if(!"organoids_num.of.objects" %in% features_to_keep) 
  features_to_keep = c(features_to_keep, "organoids_num.of.objects")
if(!"organoids_x.0.s.area_expected" %in% features_to_keep) 
  features_to_keep = c(features_to_keep, "organoids_x.0.s.area_expected")
if(!"Total.Biomass" %in% features_to_keep) 
  features_to_keep = c(features_to_keep, "Total.Biomass")

features = features[,features_to_keep]
features_metadata$CELL.LINE = substr(rownames(features_metadata), 1, 7)
```

# Stability Selection

The features selected via between-replicate correlation may still be highly redundant and unnecessary. To further reduce the number of necessary features, I perform a stepwise feature selection over the entire screen, as per Fischer et al. 2015 ("A Map of Directional Genetic Interactions in a Metazoan Cell") and Laufer et al. 2013 ("Mapping genetic interactions in human cancer cells"). This entails the following steps:
1. Choose an initial feature set. This will be the number of organoids and the total biomass, both of which suitably separated the positive from negative controls. This feature set is chosen for its biological iterpretability.
2. Model each replicate of the remaining features as a function of the already chosen feature set. Since I am interested in extracting features with the highest signal-to-noise ratio, the feature with the highest correlation between replicate residuals is added to the feature set.
3. Repeat step 2 with the new feature set until the distribution of the correlation of all remaining residuals is centered around 0, i.e. none of the remaining features add any information.

## Over Entire Screen


```r
# If data already exists, skip this part.
if(!file.exists("QC_mouse__stability_selection__features.rds")){
  source("../functions/StabilitySelection.R")
  
  # To reduce the runtime of the feature selection, I remove strongly correlating features (r >= +/-0.99).
  features_cor = cor(features, use = "pairwise.complete.obs")
  feature_cor_vector = features_cor[upper.tri(features_cor, diag = FALSE)]
  features_cor[upper.tri(features_cor, diag = TRUE)] = 0
  features_to_keep = names(which(!apply(features_cor,2,function(x) any(abs(x) > 0.99))))
  if(!"organoids_num.of.objects" %in% features_to_keep)
    features_to_keep = c(features_to_keep, "organoids_num.of.objects")
  if(!"Total.Biomass" %in% features_to_keep)
    features_to_keep = c(features_to_keep, "Total.Biomass")
  if(!"organoids_x.0.s.area_expected" %in% features_to_keep)
    features_to_keep = c(features_to_keep, "organoids_x.0.s.area_expected")
  features_reduced <- features[,features_to_keep]
  
  selection = selectByStability(
    features = features_reduced, metadata = features_metadata, 
    preselect = c("organoids_num.of.objects", "Total.Biomass", "organoids_x.0.s.area_expected"), 
    Rdim = ncol(features_reduced)-3, verbose = FALSE)
  
  # Keep only the selected features until the ratio drops to 0.5
  selected_features = na.omit(selection$selected[selection$ratioPositive >= 0.5])
  
  saveRDS(object = selection, 
          file = "QC_mouse__stability_selection__results.rds")
  saveRDS(object = features[,selected_features],
          file = "QC_mouse__stability_selection__features.rds")
  saveRDS(object = features_metadata,
          file = "QC_mouse__stability_selection__features_metadata.rds")
}

ratio_positive = na.omit(readRDS("QC_mouse__stability_selection__results.rds")$ratioPositive)
ratio_color = ifelse(ratio_positive >= 0.5, "1", "2")
plot(ggplot(data = data.frame(
  "Iteration" = seq_along(ratio_positive),
  "Ratio" = ratio_positive-0.5)) +
  geom_bar(aes(x = Iteration, y = Ratio, fill = ratio_color), stat = "identity") +
  ggtitle(label = "Ratio of features with positive correlation after each iteration") + 
  theme(legend.position = "none") + scale_y_continuous(labels = function(x) x+0.5))
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

## Individually for each Cell Line

```r
source("../functions/StabilitySelection.R")
for(cell_line in cell_lines) {
  # only run this if the data doesn't exist yet
  if(file.exists(sprintf("QC_mouse__stability_selection__features_%s.rds", cell_line))) next
  
  features_reduced = features[features_metadata$CELL.LINE == cell_line,]
  metadata_reduced = features_metadata[features_metadata$CELL.LINE == cell_line,]
  
  # To reduce the runtime of the feature selection, I remove strongly correlating features (r >= +/-0.95).
  features_cor = cor(features_reduced, use = "pairwise.complete.obs")
  feature_cor_vector = features_cor[upper.tri(features_cor, diag = FALSE)]
  features_cor[upper.tri(features_cor, diag = TRUE)] = 0
  features_to_keep = names(which(!apply(features_cor,2,function(x) any(abs(x) > 0.99))))
  if(!"organoids_num.of.objects" %in% features_to_keep)
    features_to_keep = c(features_to_keep, "organoids_num.of.objects")
  if(!"Total.Biomass" %in% features_to_keep)
    features_to_keep = c(features_to_keep, "Total.Biomass")
  if(!"organoids_x.0.s.area_expected" %in% features_to_keep)
    features_to_keep = c(features_to_keep, "organoids_x.0.s.area_expected")
  features_reduced <- features_reduced[,features_to_keep]

  selection = selectByStability(
    features_reduced, metadata_reduced, 
    preselect = c("organoids_num.of.objects", "Total.Biomass", "organoids_x.0.s.area_expected"), 
    Rdim = ncol(features_reduced)-3, verbose = FALSE)
  
  # Keep only the selected features until the ratio drops to 0.5
  selected_features = na.omit(selection$selected[selection$ratioPositive >= 0.5])
  
  saveRDS(object = selection, 
          file = sprintf("QC_mouse__stability_selection__results_%s.rds", cell_line))
  saveRDS(object = features[,selected_features],
          file = sprintf("QC_mouse__stability_selection__features_%s.rds", cell_line))
  saveRDS(object = features_metadata,
          file = sprintf("QC_mouse__stability_selection__features_metadata_%s.rds", cell_line))
}
```


```r
for(cell_line in cell_lines) {
  ratio_positive = na.omit(readRDS(sprintf("QC_mouse__stability_selection__results_%s.rds", cell_line))$ratioPositive)
  ratio_color = ifelse(ratio_positive >= 0.5, "1", "2")
  plot(ggplot(data = data.frame(
    "Iteration" = seq_along(ratio_positive),
    "Ratio" = ratio_positive-0.5)) +
    geom_bar(aes(x = Iteration, y = Ratio, fill = ratio_color), stat = "identity") +
    ggtitle(
      label = "Ratio of features with positive correlation after each iteration",
      subtitle = cell_line) +
    theme(legend.position = "none") + scale_y_continuous(labels = function(x) x+0.5))
}
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-2.png)![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-3.png)![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-4.png)

The overlap of selected features is

```r
selected_features = list()
for(cell_line in cell_lines) {
  selected_features[[cell_line]] = colnames(readRDS(
    sprintf("QC_mouse__stability_selection__features_%s.rds", cell_line)))
}
selected_features_freq = as.numeric(sort(table(unlist(selected_features)), decreasing = TRUE))
selected_features_freq = data.frame(
  "Features" = seq_along(selected_features_freq), 
  "Frequency" = selected_features_freq)
ggplot(data = selected_features_freq) + geom_point(aes(x = Features, y = Frequency)) + 
  ggtitle("Frequency of Selected Features Across Cell Lines")
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)
