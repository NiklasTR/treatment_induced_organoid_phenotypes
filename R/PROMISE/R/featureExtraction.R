#' @title Extract organoid features
#' 
#' @description Loads a projection and the corresponding segmentation map and 
#' calculates the haralick features using EBImage. The organoids are not 
#' differentiated individually
#' 
#' @param plateIndir This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the features were correctly calculated and saved
#' 
#' @author Jan Sauer
#' 
#' @examples print(extractOrganoidFeatures_texture)
#' @export
extractOrganoidFeatures_texture <- function(plateIndir, platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  # Load segmap
  segmap_fn = segmentationHdf5Filename(
    filedir = file.path(segmentationdir, platename), 
    platename = platename, row = row, col = col)
  segmap_mask = h5read(segmap_fn, "mask")
  
  # Load projection
  proj_fn = projectionHdf5Filename(
    filedir = file.path(hdf5projection, platename), 
    platename = platename, row = row, col = col)
  proj = h5read(proj_fn, "images")
  
  # Calculate features for every field
  num_fields = dim(proj)[4]
  features = list()
  for(fld in seq_len(num_fields)) {
    img = Image(data = proj[,,,fld], colormode = "color")
    mask = matrix(
      as.integer(segmap_mask[,,fld] > 0), 
      nrow = nrow(segmap_mask), 
      ncol = ncol(segmap_mask))
    features[[fld]] = computeFeatures.haralick(
      x = mask, ref = normalize(img), 
      haralick.scales = c(1, 2, 4, 8, 16, 32, 64, 128, 256))
  }
  features = do.call(what = rbind, args = features)
  feature_names = colnames(features)
  
  # Write features to file for well
  feature_dir = file.path(featuresdir, platename)
  dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
  h5FileName = featuresHdf5Filename(
    filedir = feature_dir, platename = platename, 
    row = row, col = col, configdir = configdir)
  file_made = h5createFile(h5FileName)
  if(!file_made) {H5close(); return(FALSE)}
  
  dataset_made = h5createDataset(file = h5FileName, dataset = "features_texture", dims = dim(features), storage.mode = "double",
                                 chunk = c(1, ncol(features)), level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = feature_names, file = h5FileName, name = "feature_names", level = hdf5_compression)
  h5write(obj = features, file = h5FileName, name = "features_texture")
  H5close()
  return(TRUE)
}

#' @title Extract organoid features
#' 
#' @description Loads a projection and the corresponding segmentation map and 
#' calculates the full set of EBImage features as well as the texture features 
#' of only the primitive foreground / background segmentation
#' 
#' @param plateIndir This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the features were correctly calculated and saved
#' 
#' @author Jan Sauer
#' 
#' @examples print(extractOrganoidFeatures)
#' @export
extractOrganoidFeatures <- function(plateIndir, platename, row, col, configdir) {
  # This is a necessary manual import as there seems to be a bug with EBImage:
  # computeFeatures.moment isn't loaded into the namespace. Or something. Who 
  # knows. R is a terrible language.
  library(EBImage)
  source(file.path(configdir, "watchdogConfig.R"))
  # Load segmap
  segmap_fn = segmentationHdf5Filename(
    filedir = file.path(segmentationdir, platename), 
    platename = platename, row = row, col = col)
  segmap_mask = h5read(segmap_fn, "mask")
  
  # Load projection
  proj_fn = projectionHdf5Filename(
    filedir = file.path(hdf5projection, platename), 
    platename = platename, row = row, col = col)
  proj = h5read(proj_fn, "images")
  
  # Calculate features for every field
  num_fields = dim(proj)[4]
  features = list()
  features_haralick = list()
  for(fld in seq_len(num_fields)) {
    img = Image(data = proj[,,,fld], colormode = "color")
    mask = segmap_mask[,,fld]
    mask_haralick = matrix(
      as.integer(segmap_mask[,,fld] > 0), 
      nrow = nrow(segmap_mask), 
      ncol = ncol(segmap_mask))
    
    edges = mask == 1
    edges = matrix(as.integer(edges), nrow = nrow(mask), ncol = ncol(mask))
    edges = closing(x = edges, kern = makeBrush(size = 5, shape = "disc"))
    
    # Remove open edge segments without any insides
    # Because 'fillHull' doesn't work well with image boundaries, detecting 
    # proper 'insides' of organoids requires a few extra steps
    labels = bwlabel(1 - edges)
    # The background most likely corresponds to the label with the largest volume
    labels_freq = table(as.numeric(labels))
    bg_val = as.integer(names(which.max(labels_freq)))
    fg = (labels != bg_val) * (labels != 0)
    fg = dilate(x = fg, kern = makeBrush(size = 15, shape = "disc"))
    fg = fillHull(fg)
    fg_labels = bwlabel(fg)
    # Propagate the foreground into the entire mask and set all pixels to 0 
    # which do not correspond to an organoid
    prop = propagate(
      x = matrix(0, nrow=2048, ncol=2048), seeds = fg_labels, mask = mask > 0)
    mask[prop == 0] = 0
    # Because the segmentation isn't perfect, additional watershedding helps separate
    # connected organoids in the labelling. The blurring allows for a gradient-based
    # watershedding
    boolmask = mask
    boolmask[mask > 1] = 1
    boolmask = fillHull(boolmask)
    f = makeBrush(size=151, shape = "disc")
    f = f/sum(f)
    density_mask = filter2(x = boolmask, filter = f)
    density_mask[boolmask == 0] = 0
    labeled_mask = watershed(density_mask, tolerance = 0.1)
    
    features[[fld]] = computeFeatures(
      x = labeled_mask, ref = normalize(img), 
      haralick.scales = c(1, 2, 4, 8, 16, 32, 64, 128, 256))
    
    features_haralick[[fld]] = computeFeatures.haralick(
      x = mask_haralick, ref = normalize(img), 
      haralick.scales = c(1, 2, 4, 8, 16, 32, 64, 128, 256))
  }
  features = do.call(what = rbind, args = features)
  features_haralick = do.call(what = rbind, args = features_haralick)
  feature_names = colnames(features)
  feature_names_haralick = colnames(features_haralick)
  
  # Write features to file for well
  feature_dir = file.path(featuresdir, platename)
  dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
  h5FileName = featuresHdf5Filename(
    filedir = feature_dir, platename = platename, 
    row = row, col = col, configdir = configdir)
  file_made = h5createFile(h5FileName)
  if(!file_made) {H5close(); return(FALSE)}
  
  dataset_made = h5createDataset(
    file = h5FileName, dataset = "features", dims = dim(features), 
    storage.mode = "double", chunk = c(1, ncol(features_haralick)), 
    level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  dataset_made = h5createDataset(
    file = h5FileName, dataset = "features_texture", dims = dim(features_haralick), 
    storage.mode = "double", chunk = c(1, ncol(features_haralick)), 
    level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = feature_names, file = h5FileName, name = "feature_names", level = hdf5_compression)
  h5write(obj = features, file = h5FileName, name = "features")
  h5write(obj = features_haralick, file = h5FileName, name = "features_texture")
  H5close()
  return(TRUE)
}

#' @title Combine well features
#' 
#' @description Combines the features for each individual well of a plate. 
#' This will only combine features if the features for every well have been 
#' calculated. The output is written to a file.
#' 
#' @param platename The plate that the well is on
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the features were correctly combined
#' 
#' @author Jan Sauer
#' 
#' @examples print(combineWellFeatures)
#' @export
combineWellFeatures <- function(platename, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  # Make sure all files are present
  feature_files = sort(list.files(file.path(featuresdir, platename)))
  all_wells = wells(nrWells = nrWells)
  expected_feature_files = character(nrWells)
  for(i in seq_len(nrow(all_wells))) {
    r = all_wells[i, "rows"]
    c = all_wells[i, "cols"]
    expected_feature_files[i] = featuresHdf5Filename(
      filedir = file.path(featuresdir, platename), platename = platename, 
      row = r, col = c, configdir = configdir)
  }
  if(!identical(feature_files, basename(expected_feature_files))) {
    warning(sprintf("Not all feature files present for '%s'", platename))
    return(FALSE)
  }
  
  features = setNames(
    object = vector(mode = "list", length = nrWells), 
    nm = all_wells$names)
  feature_names = setNames(
    object = vector(mode = "list", length = nrWells), 
    nm = all_wells$names)
  well_names = c()
  for(i in seq_len(nrow(all_wells))) {
    r = all_wells[i, "rows"]
    c = all_wells[i, "cols"]
    n = all_wells[i, "names"]
    full_path = featuresHdf5Filename(
      filedir = file.path(featuresdir, platename), platename = platename, 
      row = r, col = c, configdir = configdir)
    features[[n]] = data.frame(h5read(file = full_path, name = "features"))
    feature_names[[n]] = h5read(file = full_path, name = "feature_names")
    colnames(features[[n]]) = feature_names[[n]]
    well_names = c(well_names, rep_len(n, nrow(features[[n]])))
    H5close()
  }
  
  # Test that all feature names are identical
  if(length(unique(feature_names)) != 1) {
    warning(sprintf(
      "Feature names for '%s' are not identical across all wells", platename))
    return(NULL)
  }
  
  # Combine features
  features = do.call(rbind, features)
  
  # Sanity check
  if(length(well_names) != nrow(features)) {
    warning(sprintf(
      "Something went wrong when creating the well names vector for '%s'", 
      platename))
    return(NULL)
  }
  
  # Move all existing well features to a subfolder
  dir.create(path = file.path(
    featuresdir, platename, "wells"), recursive = FALSE, showWarnings = FALSE)
  for(fn in feature_files) {
    from = file.path(featuresdir, platename, fn)
    to = file.path(featuresdir, platename, "wells", fn)
    file.rename(from = from, to = to)
  }
  
  # Save plate features
  h5FileName = file.path(
    featuresdir, platename, sprintf("%s_features.h5", platename))
  file_made = h5createFile(h5FileName)
  if(!file_made) {H5close(); return(FALSE)}
  dataset_made = h5createDataset(
    file = h5FileName, dataset = "features", dims = dim(features), 
    storage.mode = "double", chunk = c(1, ncol(features)), 
    level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = colnames(features), file = h5FileName, 
          name = "feature_names", level = hdf5_compression)
  h5write(obj = well_names, file = h5FileName, 
          name = "well_names", level = hdf5_compression)
  h5write(obj = as.matrix(features), file = h5FileName, name = "features")
  H5close()
  return(TRUE)
}

#' @title Clean up segmentation mask
#' 
#' @description "Cleans up" a segmentation mask. This means that two separate 
#' masks will be created: one mask containing all closed organoids and one 
#' mask containing all foreground and edge fragments not inside these closed 
#' organoids. Features for these two classes of items are calculated 
#' separately.
#' 
#' @param segmap A single 2D matrix with integer values
#' 
#' @return A named list of 2D matrices. Currently with the entries "organoids" 
#' and "other".
#' 
#' @author Jan Sauer
#' 
#' @examples print(cleanupSegmap)
#' @export
cleanupSegmap <- function(segmap) {
  warning("This function is not yet fully implemented!")
  return(NULL)
  
  if(length(dim(segmap)) != 2) {
    warning("'segmap' must be a 2D matrix-like")
    return(NULL)
  }
  edges = segmap == 1
  edges = matrix(as.integer(edges), nrow = nrow(segmap), ncol = ncol(segmap))
  edges = closing(x = edges, kern = makeBrush(size = 5, shape = "disc"))
  
  # Because 'fillHull' doesn't work well with image boundaries, detecting 
  # proper 'insides' of organoids requires a few extra steps
  # 1. Isolate the background
  
  labels = bwlabel(1 - edges)
  # The background most likely corresponds to the label with the largest volume
  labels_freq = table(as.numeric(labels))
  bg_val = as.integer(names(which.max(labels_freq)))
  fg = (labels != bg_val) * (labels != 0)
  fg = closing(x = fg, kern = makeBrush(size = 5, shape = "disc"))
  # edgesfg = edges + fg
  # edgesfg[edgesfg > 1] = 1
  fg = fillHull(fg)
  organoids = matrix(0, nrow = nrow(segmap), ncol = ncol(segmap))
  organoids[edges == 1] = 1
  organoids[fg == 1] = 2
}

#' @title Calculates the mean feature vector for an entire well
#'
#' @description Calculates the mean feature vector as well as the standard deviation of all features over all organoids for a single well.
#' Also determines the mean number of cells in each organoid as well as how many organoids and single cells there are.
#'
#' @param plateIndir This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#'
#' @return The combined feature vector of the well
#'
#' @author Jan Sauer
#'
#' @examples print(calcFeaturesPerWell)
#' @export
calcFeaturesPerWell <- function(plateIndir, platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  features_filename = featuresHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = row, col = col, configdir = configdir)
  mask_filename = featuresMaskHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = row, col = col, configdir = configdir)
  
  # Read features and masks
  features = h5read(features_filename, name = "features")
  feature_names = h5read(features_filename, name = "feature_names")
  metadata = h5read(features_filename, name = "metadata")
  H5close()
  masks = h5read(mask_filename, name = "masks")
  H5close()

  # Calculate mean and standard deviation of features and compress into one vector
  features_mean = apply(features, 2, function(x) mean(x, trim = 0.2))
  features_sd = apply(features, 2, function(x) sd(x))
  well_features = c(features_mean, features_sd)
  names(well_features) = c(paste0("mean__", feature_names), paste0("sd__", feature_names))
  
  # Determine how many organoids are in the well and how many single cells outside of organoids are present
  num_organoids = nrow(features)
  single_cells = (1 - masks[,,1,]) * masks[,,3,]
  single_cells_label = bwlabel(single_cells)
  num_cells = sum(apply(single_cells_label, 3, max))
  well_features = c(well_features, Num_Organoids = num_organoids, Num_Isolated_Cells = num_cells)
  
  # DIAGNOSTIC FEATURE: Sometimes, there is a FITC Signal without a DAPI Signal. This probably shouldn't be the case. 
  # This feature looks at isolated FITC Signals, i.e. without corresponding DAPI or Cy3 signals
  # For this, the Cy3 and DAPI blots are grown to account for dye bleeding
  cy3_mask_grow = dilate(masks[,,1,], kern = makeBrush(size = 25, shape = "diamond"))
  dapi_mask_grow = dilate(masks[,,3,], kern = makeBrush(size = 25, shape = "diamond"))
  isolated_fitc = (1 - cy3_mask_grow) * (1 - dapi_mask_grow) * masks[,,2,]
  isolated_fitc_label = bwlabel(isolated_fitc)
  num_isolated_fitc = sum(apply(isolated_fitc_label, 3, max))
  size_isolated_fitc = sum(isolated_fitc != 0)
  well_features = c(well_features, Num_Isolated_FITC_Signals = num_isolated_fitc, Area_Isolated_FITC_Signals = size_isolated_fitc)
  
  return(well_features)
}

#' @title Calculates the mean feature vector for an entire plate
#'
#' @description Runs calcFeaturesPerWell() for each well and saves it as a file for the plate
#'
#' @param plateIndir This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param platename The plate that the well is on
#' @param row This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param col This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param configdir The directory of the configuration files
#'
#' @return A boolean value indicating if the function was successful
#'
#' @author Jan Sauer
#'
#' @examples print(calcFeaturesPerPlate)
#' @export
calcFeaturesPerPlate <- function(plateIndir, platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))

  W = wells(nrWells)

  fn = featuresHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = W$rows[1], col = W$cols[1], configdir = configdir)
  metadata = h5read(file = fn, name = "metadata")
  metadata = metadata[which(metadata[,1] != "Well"),]
  metadata[metadata[,1] == "PackageVersion",2] = as.character(packageVersion("PROMISE"))
  
  features = list()
  for(i in seq_len(nrow(W))) {
    cat(i, "\n")
    features[[W[i,"names"]]] = calcFeaturesPerWell(platename = platename, row = W[i,"rows"], col = W[i,"cols"], configdir = configdir)
  }
  
  features = do.call(rbind, features)
  
  h5filename = file.path(featuresdir, platename, sprintf("%s_features.h5", platename))
  file_created = h5createFile(h5filename)
  if(!file_created) {H5close(); return(FALSE)}
  
  dataset_made = h5createDataset(file = h5filename, dataset = "features", dims = dim(features), storage.mode = storage.mode(features),
                                 chunk = c(1, ncol(features)), level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  
  h5write(obj = features, file = h5filename, name = "features")
  h5write(obj = rownames(features), file = h5filename, name = "wells", level = hdf5_compression)
  h5write(obj = colnames(features), file = h5filename, name = "feature_names", level = hdf5_compression)
  h5write(obj = metadata, file = h5filename, name = "metadata", level = hdf5_compression)
  H5close()
  return(TRUE)
}

#' @title Analyses the per-plate features
#'
#' @description Performs quality control on the per-plate features followed by a a shrinkage linear discriminant analysis (bioconductor package 'sda').
#' If there are multiple plates, use these as replicates. If there is only one plate, however, use wells with identical treatment on the plate as
#' biological replicates. The current implementation requires that there is an equal number of wells for each drug treatment
#'
#' @param platenames A list of platenames to use as biological replicates for the analysis. This list can be one element long
#' @param configdir The directory of the configuration files
#' @param grouping Can take the values "drugs", "combo", or "type". "Drugs" indicates that the drugs alone should be used to label the wells, 
#' "combo" indicates that the combination of drugs and concentration should be used, and "type" indicates that the treatment type should be used
#' ("drug", "positive control", or "negative control")
#' @param correlation_threshold Only features, which correlate between replicates are used for the analysis. This threshold indicates how strongly
#' they should correlate
#'
#' @return A boolean value indicating if the function was successful
#'
#' @author Jan Sauer
#'
#' @examples print(analysePlateFeatures)
#' @export
analysePlateFeatures <- function(platenames, configdir, grouping, correlation_threshold) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  # Make sure the plates are actual replicates
  if(length(unique(substr(platenames, 1, 7))) != 1) stop("Biological replicates must belong to the same patient and sample")
  if(length(unique(substr(platenames, 12, 14))) != 1) stop("Biological replicates must have the same layout")
  if(!grouping %in% c("drugs", "combo", "type")) stop("'grouping' must be either 'drugs', 'combo', or 'type'")
  
  # Load plate layout
  layoutID = unique(substr(platenames, 12, 14))
  layout = read.table(file.path(layoutdir, sprintf("%s.txt", layoutID)), header = TRUE, stringsAsFactors = FALSE)
  if(grouping == "drugs") layout$treatment = layout$substance
  if(grouping == "combo") layout$treatment = paste0(layout$substance, "_", layout$concentration)
  if(grouping == "type") layout$treatment = layout$group
  
  # Load features
  features = list()
  for(platename in platenames) {
    features_fn = file.path(featuresdir, platename, sprintf("%s_features.h5", platename))
    features[[platename]] = h5read(file = features_fn, name = "features")
    feature_names = h5read(file = features_fn, name = "feature_names")
    wells = h5read(file = features_fn, name = "wells")
    metadata = h5read(file = features_fn, name = "metadata")
    rownames(features[[platename]]) = wells
    colnames(features[[platename]]) = feature_names
  
    # Match the order of wells in the layout with the order of wells in the features
    features[[platename]] = features[[platename]][match(layout$well, rownames(features[[platename]])),]
  }

  # Quality control
  # Check that the distribution of values for each feature is "continuous", i.e. at least 1/3 of the wells have unique values for the features
  for(platename in platenames) {
    unique_vals = apply(features[[platename]], 2, function(x) length(unique(x)))
    features[[platename]] = features[[platename]][, unique_vals >= nrow(features[[platename]]) / 3]
  }
  
  # Use only the features present in all replicates (if replicates are present)
  if(length(platenames) > 1) {
    all_feature_names = lapply(features, colnames)
    common_feature_names = do.call(intersect, setNames(all_feature_names, NULL))
    features = lapply(features, function(x) x[,common_feature_names])
  }
  
  # CORRELATION
  # Filter out features not correlated between replicates
  # If there are multiple plates, calculate the correlation of each feature between plates across all wells.
  # If there is only one plate, calculate the correlation across wells with identical treatments. For this, treat
  # each plate as a combination of replicates and generate N subsets of wells representing the replicates
  # Calculate the mean correlation between all replicates (if more than 2).
  if(length(platenames) > 1) {
    # Find the mean values over the treatments (depending on 'grouping')
    features_treatments = lapply(features, function(x) {
      mat = aggregate(x, list(G = layout$treatment), mean)
      rownames(mat) = mat$G
      mat[,-1]
    })
    correlations = list()
    for(i in seq_len(length(features)-1)) {
      for(j in seq_along(features)[-seq_len(i)]) {
        correlations = c(correlations, list(diag(cor(features_treatments[[i]], features_treatments[[j]]))))
      }
    }
    correlations = apply(do.call(rbind, correlations), 2, mean)
  } else {
    num_replicates = unique(table(layout[layout$group == "drug", "treatment"]))
    treatments = unique(layout[layout$group == "drug", "treatment"])
    if(length(num_replicates) != 1) {
      stop("Correlation across replicates cannot be calculated because there is an uneven distribution of on-plate replicates")
    }
    replicates = matrix(nrow = 0, ncol = num_replicates)
    for(treatment in treatments) {
      rep_wells = layout[layout$treatment == treatment, "well"]
      replicates = rbind(replicates, rep_wells)
    }
    rep_features = list()
    for(rep in seq_len(num_replicates)) {
      rep_features[[rep]] = features[[platename]][replicates[,rep],]
    }
    
    correlations = list()
    for(i in seq_len(num_replicates-1)) {
      for(j in seq_len(num_replicates)[-seq_len(i)]) {
        correlations = c(correlations, list(diag(cor(rep_features[[i]], rep_features[[j]]))))
      }
    }
    correlations = apply(do.call(rbind, correlations), 2, mean)
  }
  
  # Use only the features with a replicate correlation of >= correlation_threshold
  features = lapply(features, function(x) x[,names(which(correlations >= correlation_threshold))])
  
  # Combine the replicates into a single matrix
  features = do.call(rbind, lapply(features, data.frame))
  labels_features = rep(layout$treatment, times=length(platenames))
  type_features = rep(layout$group, times=length(platenames))
  
  # # SDA
  # labels = layout[match(layout$well, rownames(features)),"treatment"]
  # model = sda(features, labels, lambda = 1)
  # psda = predict.sda(model, features)
  
  # LDA
  # # Train the classifier on a training set and validate with a test set
  # training_id = sample(x = seq_len(nrow(features)), size = 3 * nrow(features) / 4, replace = FALSE)
  # testing_id = seq_len(nrow(features))[-training_id]
  # model.lda = MASS::lda(features[training_id,], labels_features[training_id])
  # plda = predict(model.lda, features[testing_id,])
  # accuracy = sum(plda$class == labels_features[testing_id]) / length(testing_id)
  
  # Train the classifier on the entire data set to look at clustering
  # LDA requires that features are not constant within groups.
  split_features = split(features, layout$treatment)
  split_sd = lapply(split_features, function(x) apply(x, 2, sd))
  split_sd = do.call(rbind, split_sd)
  min_split_sd = apply(split_sd, 2, min)
  features = features[,min_split_sd > 0]
  
  model.lda = MASS::lda(features, labels_features)
  var_explained = round(model.lda$svd^2 / sum(model.lda$svd^2)*100, digits = 2)
  plda = predict(model.lda, features)
  plotdat = data.frame(class = labels_features, type = type_features, lda = plda$x)
  gpl = ggplot(plotdat) + geom_point(aes(lda.LD1, lda.LD2, colour = class, shape = type), size = 2.5) +
    labs(x = paste("LD1 (", var_explained[1], "%)", sep=""),
         y = paste("LD2 (", var_explained[2], "%)", sep=""))
  
  # Calculate variances of plotdat
  # within-groups
  variances = c()
  for(group in unique(plotdat$class)) {
    subdat = as.matrix(plotdat[plotdat$class == group, c(-1,-2)])
    submean = colMeans(subdat)
    meansmat = matrix(rep(submean, nrow(subdat)), nrow = nrow(subdat), byrow = TRUE)
    variances = c(variances, sum(sqrt(rowSums((subdat - meansmat)^2))))
  }
  withinvar = sum(variances) / nrow(plotdat)
  
  # total
  totmean = colMeans(plotdat[,-c(1,2)])
  totmeansmat = matrix(rep(totmean, nrow(plotdat)), nrow = nrow(plotdat), byrow = TRUE)
  totvar = mean(sqrt(rowSums((plotdat[,-c(1,2)] - totmeansmat)^2)))
  
  return(list(data = plotdat, plot = gpl, withinvar = withinvar, betweenvar = totvar - withinvar, totalvar = totvar))
}

#' @title Combines the features per well into a single file for an entire plate
#' 
#' @description Combines all features for each individual well of a plate into a single file. This function requires that 
#' features for well A01 exist for the extraction of feature names.
#' 
#' @param platename The plate that the well is on
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the combining was successful
#' 
#' @author Jan Sauer
#' 
#' @examples print(combineFeatures)
#' @export
combineFeatures <- function(platename, configdir) {
  warning("DEBUG: 'combineFeatures()' is currently not maintained. Use with caution! It's still in the package because it may be useful later.")
  return(NULL)
  source(file.path(configdir, "watchdogConfig.R"))
  
  inputFilename = featuresHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = "A", col = "01", configdir = configdir)
  features_names = h5read(inputFilename, name = "feature_names")
  H5close()
  
  features_list = list()
  wells_vector = c()
  W = wells(nrWells)
  for(i in seq_len(nrow(W))) {
    r = W[i, "rows"]
    c = W[i, "cols"]
    well = W[i, "names"]
    inputFilename = featuresHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = r, col = c, configdir = configdir)
    dat = h5read(inputFilename, name = "features")
    features_list[[i]] = dat
    wells_vector = c(wells_vector, rep(well, nrow(dat)))
    H5close()
  }
  features = do.call(rbind, features_list)
  
  h5filename = file.path(featuresdir, platename, sprintf("%s_features.h5", platename))
  file_created = h5createFile(h5filename)
  if(!file_created) {H5close(); return(FALSE)}
  
  dataset_made = h5createDataset(file = h5filename, dataset = "features", dims = dim(features), storage.mode = "double",
                                 chunk = c(1, ncol(features)), level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  dataset_made = h5createDataset(file = h5filename, dataset = "wells", dims = length(wells_vector), storage.mode = "character", 
                                 size = 10, level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  
  h5write(features, h5filename, "features")
  h5write(wells_vector, h5filename, "wells")
  H5close()
  return(TRUE)
}

#' @title Extract organoid features
#' 
#' @description !THIS METHOD IS OBSOLETE! Segments the projected images via thresholding to obtain a rough estimate of the organoids and then calculates the 
#' features (EBImage::computeFeatures) for each individual organoid. Both the mask and the features are saved as an hdf5 file for each well.
#' 
#' @param plateIndir This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the features were correctly calculated and saved
#' 
#' @author Jan Sauer
#' 
#' @examples print(extractOrganoidFeatures_v1)
#' @export
extractOrganoidFeatures_v1 <- function(plateIndir, platename, row, col, configdir) {
  warning("This version of the 'extractOrganoidFeatures()' uses a thresholding-based method of segmentation. The newer implementation combines this with a sobel edge detector.")
  return(NULL)
  
  source(file.path(configdir, "watchdogConfig.R"))
  
  # This is the overlap of the microscope images. Any organoid closer than overlap/2 pixels to an edge with a neighboring field is discarded (as it'll be counted in the neighboring field)
  overlap = 330
  
  # Load input images
  inputFile = projectionHdf5Filename(filedir = file.path(hdf5projection, platename), platename = platename, row = row, col = col, configdir = configdir)
  img = h5read(inputFile, name="images")
  metadata = h5read(inputFile, name="metadata")
  metadata[metadata[,1] == "PackageVersion",2] = as.character(packageVersion("PROMISE"))
  
  # Smooth images
  filter_dapi = makeBrush(size=51, shape="Gaussian", sigma = 1)
  filter_cy3_fitc = makeBrush(size=51, shape="Gaussian", sigma = 3)
  
  features_list = list()
  fields_list = list()
  for(fld in seq_len(dim(img)[4])) {
    img_cy3_smooth = filter2(x = img[,,1,fld], filter = filter_cy3_fitc)
    img_fitc_smooth = filter2(x = img[,,2,fld], filter = filter_cy3_fitc)
    img_dapi_smooth = filter2(x = img[,,3,fld], filter = filter_dapi)
    
    # Find organoids
    img_bin = thresh(x = normalize(img_cy3_smooth), w = 250, h = 250, offset = 0.05)
    img_bin = fillHull(img_bin)
    img_bin = opening(img_bin, kern = makeBrush(5))
    img_label = bwlabel(img_bin)
    
    # Find nuclei
    nucleus_thresh = thresh(x = normalize(img_dapi_smooth), w = 10, h = 10, offset = 0.02)
    nucleus_opening = opening(x = nucleus_thresh)
    
    # Find dead nuclei
    dead_nucleus_thresh = thresh(x = normalize(img_fitc_smooth), w = 10, h = 10, offset = 0.02)
    dead_nucleus_opening = opening(x = dead_nucleus_thresh)
    
    # Match the number of nuclei to each organoid
    num_cells = c()
    for(i in seq_len(max(img_label))) {
      mask = img_label == i
      num_cells = c(num_cells, max(bwlabel(nucleus_opening * mask)))
    }
    
    # I've added this 'require' statement because otherwise R throws an error 'computeFeatures.moment not found'.
    require(EBImage)
    features = computeFeatures(x = img_label, ref = list(img[,,1,fld], img[,,2,fld], img[,,3,fld]), haralick.scales=c(1, 2, 4, 8, 16, 32, 64, 128), refnames = c("Cy3", "FITC", "DAPI"))
    features = cbind(features, "Num_Cells_Per_Organoid" = num_cells)
    
    # Discard organoids based on edge proximity
    indices_to_remove = c()
    for(i in seq_len(nrow(features))) {
      x = features[i,"x.0.m.cx"]
      y = features[i,"x.0.m.cy"]
      if(x < overlap/2 | x > (2048 - overlap/2) | y < overlap/2 | y > (2048 - overlap/2)) indices_to_remove = c(indices_to_remove, i)
    }
    if(length(indices_to_remove) > 0) features = features[-indices_to_remove,]
    
    features_list[[fld]] = features
    fields_list[[fld]] = combine(img_bin, dead_nucleus_opening, nucleus_opening)
  }
  
  features = do.call(what = "rbind", args = features_list)
  rownames(features) = NULL
  
  fields_array = abind(fields_list, along = 4)
  
  # Save image masks as hdf5 file
  dir.create(file.path(featuresdir, platename), showWarnings = FALSE, recursive = TRUE)
  h5maskFileName = featuresMaskHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = row, col = col, configdir = configdir)
  file_made = h5createFile(h5maskFileName)
  if(!file_made) {H5close(); return(FALSE)}
  dataset_made = h5createDataset(file = h5maskFileName, dataset = "masks", dims = dim(fields_array), H5type = "H5T_NATIVE_UINT16",
                                 chunk = c(hdf5_chunksize, hdf5_chunksize, 1, 1), level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = metadata, file = h5maskFileName, name = "metadata", level = hdf5_compression)
  h5write(obj = fields_array, file = h5maskFileName, name = "masks")
  H5close()
  
  # Create hdf5 file (and directory if necessary)
  h5FileName = featuresHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = row, col = col, configdir = configdir)
  file_made = h5createFile(h5FileName)
  if(!file_made) {H5close(); return(FALSE)}
  
  dataset_made = h5createDataset(file = h5FileName, dataset = "features", dims = dim(features), storage.mode = "double",
                                 chunk = c(1, ncol(features)), level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = metadata, file = h5FileName, name = "metadata", level = hdf5_compression)
  h5write(obj = colnames(features), file = h5FileName, name = "feature_names", level = hdf5_compression)
  h5write(obj = features, file = h5FileName, name = "features")
  H5close()
  return(TRUE)
}

#' @title Extract organoid features
#' 
#' @description !THIS METHOD IS OBSOLETE! Segments the projected images via thresholding and edge detection to obtain a rough estimate of the organoids and then calculates the 
#' features (EBImage::computeFeatures) for each individual organoid. Both the mask and the features are saved as an hdf5 file for each well.
#' 
#' @param plateIndir This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the features were correctly calculated and saved
#' 
#' @author Jan Sauer
#' 
#' @examples print(extractOrganoidFeatures_v2)
#' @export
extractOrganoidFeatures_v2 <- function(plateIndir, platename, row, col, configdir) {
  warning("This version of 'extractOrganoidFeatures()' has primitive thresholding. The newer version uses the DNN segmap")
  return(NULL)
  source(file.path(configdir, "watchdogConfig.R"))
  
  # This is the overlap of the microscope images. Any organoid closer than overlap/2 pixels to an edge with a neighboring field is discarded (as it'll be counted in the neighboring field)
  overlap = 105
  
  # Load input images
  inputFile = projectionHdf5Filename(filedir = file.path(hdf5projection, platename), platename = platename, row = row, col = col, configdir = configdir)
  img = h5read(inputFile, name="images")
  metadata = h5read(inputFile, name="metadata")
  metadata[metadata[,1] == "PackageVersion",2] = as.character(packageVersion("PROMISE"))
  
  # Smoothing filters
  filter_dapi = makeBrush(size=51, shape="Gaussian", sigma = 1)
  filter_cy3_broad = makeBrush(size=51, shape="Gaussian", sigma = 9)
  filter_cy3_small = makeBrush(size=51, shape="Gaussian", sigma = 3)
  
  # Edge filters
  gradient_filters = list(
    matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3, ncol = 3),
    matrix(c(-1, -2, -1, 0, 0, 0, 1, 2, 1), nrow = 3, ncol = 3),
    matrix(c(2, 1, 0, 1, 0, -1, 0, -1, -2), nrow = 3, ncol = 3),
    matrix(c(-2, -1, 0, -1, 0, 1, 0, 1, 2), nrow = 3, ncol = 3),
    matrix(c(1, 0, -1, 2, 0, -2, 1, 0, -1), nrow = 3, ncol = 3),
    matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), nrow = 3, ncol = 3),
    matrix(c(0, 1, 2, -1, 0, 1, -2, -1, 0), nrow = 3, ncol = 3),
    matrix(c(0, -1, -2, 1, 0, -1, 2, 1, 0), nrow = 3, ncol = 3))
  
  features_list = list()
  fields_list = list()
  for(fld in seq_len(dim(img)[4])) {
    # Detect actin backbone structure using both a window thresholding and an edge detector
    # NOTE: This is too unreliable in its current implementation.
    #       The idea is to do a voronoi expansion from each backbone onto the nuclei to determine
    #       which nuclei belong to which organoid
    # img_cy3_smooth_broad = filter2(x = normalize(img[,,1,fld]), filter = filter_cy3_broad, boundary = 0)
    # img_cy3_smooth_broad = img_cy3_smooth_broad / quantile(img_cy3_smooth_broad, 0.99)
    # img_cy3_bin_backbone = thresh(x = img_cy3_smooth_broad, w = 5, h = 5, offset = 0.0025)
    # img_cy3_bin_backbone = closing(img_cy3_bin_backbone, kern = makeBrush(9, shape="disc"))
    
    # Perform edge detection on cy3 channel
    img_cy3_smooth_fine = filter2(x = normalize(img[,,1,fld]), filter = filter_cy3_small, boundary = 0)
    img_cy3_smooth_fine = img_cy3_smooth_fine / quantile(img_cy3_smooth_fine, 0.99)
    edges = matrix(0, nrow = nrow(img_cy3_smooth_fine), ncol = ncol(img[,,1,fld]))
    for(gf in gradient_filters) {
      edges = pmax(edges, filter2(normalize(img_cy3_smooth_fine), gf, boundary = "replicate"))
    }
    edges_bin = edges > 0.025
    
    # Detect the likely location of organoids, i.e. organoids than are detected by the thresholding
    img_organoids = thresh(x = normalize(img_cy3_smooth_fine), w = 250, h = 250, offset = 0.05)
    
    # Combine thresholding with edge detection
    img_bin = img_organoids + edges_bin > 0
    img_bin = closing(img_bin, kern = makeBrush(21, shape="disc"))
    
    # Fill the objects (fillHull() doesn't properly take image boundaries into consideration, so this alternative method is required)
    img_bin_inv = 1 - img_bin
    img_bin_inv_label = bwlabel(img_bin_inv)
    sizes = table(img_bin_inv_label)
    # '0' corresponds to the foreground of the original image
    # The largest area of the inverted image SHOULD correspond to the true background of the original binary image
    sizes_no0 = sizes[names(sizes) != "0"]
    true_background = names(sizes_no0)[which.max(sizes_no0)]
    # Set all other entries (not the true background) to 0 on the inverted image (equivalent to 1 on the original binary image)
    img_bin_inv_label[img_bin_inv_label != as.numeric(true_background)] = 0
    img_bin = 1 - img_bin_inv_label
    
    img_bin = opening(img_bin, kern = makeBrush(21, shape="disc"))
    
    # Detect cells
    img_dapi_smooth = filter2(x = normalize(img[,,3,fld] + img[,,2,fld]), filter = filter_dapi, boundary = 0)
    img_dapi_smooth = img_dapi_smooth / quantile(img_dapi_smooth, 0.99)
    nucleus_thresh = thresh(x = normalize(img_dapi_smooth), w = 10, h = 10, offset = 0.02)
    nucleus_opening = opening(x = nucleus_thresh)
    
    # Find dead nuclei
    img_fitc_smooth = filter2(x = img[,,2,fld], filter = filter_dapi)
    img_fitc_smooth = img_fitc_smooth / quantile(img_fitc_smooth, 0.99)
    dead_nucleus_thresh = thresh(x = normalize(img_fitc_smooth), w = 10, h = 10, offset = 0.02)
    dead_nucleus_opening = opening(x = dead_nucleus_thresh)
    
    # Match the number of nuclei to each organoid
    img_label = bwlabel(img_bin)
    num_cells = c()
    for(i in seq_len(max(img_label))) {
      mask = img_label == i
      num_cells = c(num_cells, max(bwlabel(nucleus_opening * mask)))
    }
    
    # I've added this 'require' statement because otherwise R throws an error 'computeFeatures.moment not found'.
    require(EBImage)
    features = computeFeatures(x = img_label, ref = list(img[,,1,fld], img[,,2,fld], img[,,3,fld]), haralick.scales=c(1, 2, 4, 8, 16, 32, 64, 128), refnames = c("Cy3", "FITC", "DAPI"))
    features = cbind(features, "Num_Cells_Per_Organoid" = num_cells)
    
    # Discard organoids based on edge proximity
    indices_to_remove = c()
    for(i in seq_len(nrow(features))) {
      x = features[i,"x.0.m.cx"]
      y = features[i,"x.0.m.cy"]
      if(x < overlap/2 | x > (2048 - overlap/2) | y < overlap/2 | y > (2048 - overlap/2)) indices_to_remove = c(indices_to_remove, i)
    }
    if(length(indices_to_remove) > 0) features = features[-indices_to_remove,]
    
    features_list[[fld]] = features
    fields_list[[fld]] = combine(img_bin, dead_nucleus_opening, nucleus_opening)
  }
  
  features = do.call(what = "rbind", args = features_list)
  rownames(features) = NULL
  
  fields_array = abind(fields_list, along = 4)
  
  # Save image masks as hdf5 file
  dir.create(file.path(featuresdir, platename), showWarnings = FALSE, recursive = TRUE)
  h5maskFileName = featuresMaskHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = row, col = col, configdir = configdir)
  file_made = h5createFile(h5maskFileName)
  if(!file_made) {H5close(); return(FALSE)}
  dataset_made = h5createDataset(file = h5maskFileName, dataset = "masks", dims = dim(fields_array), H5type = "H5T_NATIVE_UINT16",
                                 chunk = c(hdf5_chunksize, hdf5_chunksize, 1, 1), level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = metadata, file = h5maskFileName, name = "metadata", level = hdf5_compression)
  h5write(obj = fields_array, file = h5maskFileName, name = "masks")
  H5close()
  
  # Create hdf5 file (and directory if necessary)
  h5FileName = featuresHdf5Filename(filedir = file.path(featuresdir, platename), platename = platename, row = row, col = col, configdir = configdir)
  file_made = h5createFile(h5FileName)
  if(!file_made) {H5close(); return(FALSE)}
  
  dataset_made = h5createDataset(file = h5FileName, dataset = "features", dims = dim(features), storage.mode = "double",
                                 chunk = c(1, ncol(features)), level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = metadata, file = h5FileName, name = "metadata", level = hdf5_compression)
  h5write(obj = colnames(features), file = h5FileName, name = "feature_names", level = hdf5_compression)
  h5write(obj = features, file = h5FileName, name = "features")
  H5close()
  return(TRUE)
}
