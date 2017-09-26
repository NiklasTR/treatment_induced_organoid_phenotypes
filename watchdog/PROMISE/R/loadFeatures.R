#' @title Load all features for a plate
#' 
#' @description Retrieve all feature for a plate
#' 
#' @param platename The plate ID
#' @param configdir The directory of the configuration files
#' @param feature_type The type of feature to return. See details.
#' 
#' @details The feature_type parameter can have the following values:
#' \describe{
#'   \item{organoids}{Returns features for individual organoids. Touching 
#'   organoids are separated as best as possible.}
#'   \item{clumps}{Returns features for individual clumps. A "clump" is any 
#'   contiguous foreground patch without any object separation.} 
#'   \item{foreground}{Returns a single feature set for the entire foreground 
#'   of an image without any object separation. This feature set only includes 
#'   shape-independent haralick features}}
#' 
#' @return A data.frame of features. Each row represents one well's field, 
#'         each column represents a feature.
#' 
#' @author Jan Sauer
#' 
#' @examples print(loadFeatures)
#' @export
loadFeatures <- function(platename, configdir, feature_type) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  # hdf5 key selection
  hdf5key = switch(
    feature_type, 
    "organoids" = c("features", "feature_names", "well_names"),
    "clumps" = c("features_clumps", "feature_names_clumps", "well_names_clumps"),
    "foreground" = c("features_noseg", "feature_names_noseg", "well_names_noseg"))
  if(is.null(hdf5key)) {
    warning(
      "'feature_type' must be one of 'organoids', 'clumps', or 'foreground'")
    return(NULL)
  }
  
  feature_fn = file.path(
    featuresdir, platename, sprintf("%s_features.h5", platename))
  features = data.frame(h5read(file = feature_fn, name = hdf5key[1]))
  feature_names = as.character(h5read(file = feature_fn, name = hdf5key[2]))
  wells = as.character(h5read(file = feature_fn, name = hdf5key[3]))
  colnames(features) = feature_names
  features$WELL = wells
  
  return(features)
  
}

#' @title Load features for a well
#' 
#' @description Retrieve all features for a well
#' 
#' @param platename The plate ID
#' @param row The plate row
#' @param col The plate column
#' @param configdir The directory of the configuration files
#' 
#' @return A list of feature dataframes. The entries are as follows: 
#'   "organoid_features" correspond to features for individual organoids.
#'   "clump_features" corespond to features of entire contiguous foreground 
#'   patches without any (i.e. groups of connected organoids)
#'   "texture_features" correspond to the shape-independent haralick features 
#'   of the foreground without any object separation.
#' 
#' @author Jan Sauer
#' 
#' @examples print(loadFeaturesForWell)
#' @export
loadFeaturesForWell <- function(platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  # The feature file could be in the main folder or in the 'wells' subfolder
  # if it's already been compressed
  feature_fn = featuresHdf5Filename(
    filedir = file.path(featuresdir, platename), platename = platename, 
    row = row, col = col)
  if(!file.exists(feature_fn)) {
    feature_fn = featuresHdf5Filename(
      filedir = file.path(featuresdir, platename, "wells"), 
      platename = platename, row = row, col = col)
  }
  if(!file.exists(feature_fn)) {
    warning(sprintf(
      "Feature file for '%s / %s / %s' doesn't exist",
      platename, row, col))
    return(NULL)
  }
  
  features_organoids = data.frame(h5read(file = feature_fn, name = "features"))
  features_clumps = data.frame(h5read(file = feature_fn, name = "features_clumps"))
  features_noseg = data.frame(h5read(file = feature_fn, name = "features_noseg"))
  colnames(features_organoids) = h5read(file = feature_fn, name = "feature_names")
  colnames(features_clumps) = h5read(file = feature_fn, name = "feature_names_clumps")
  colnames(features_noseg) = h5read(file = feature_fn, name = "feature_names_noseg")
  return(
    list("organoid_features" = features_organoids, 
         "clump_features" = features_clumps, 
         "texture_features" = features_noseg))
}

#' @title Get feature vector by compound
#' 
#' @description Retrieve all feature vectors corresponding to a compound. This 
#'              searches through all plates currently found in the "featuresdir"
#'              location
#' 
#' @param compound The compound for which to retrieve the feature vectors
#' @param feature_type The subfolder of the "features" PROMISE folder from 
#'                     which to retrieve the features (e.g. "inception" for 
#'                     inception features)
#' @param organoid_type Either "mouse", "human", or "both", depending on which 
#'                      organoid data to search through
#' @param configdir The directory of the configuration files
#' 
#' @return A matrix of features. Each row represents one well, each column 
#'         represents a feature. The rownames indicate the source well
#' 
#' @author Jan Sauer
#' 
#' # @examples print(get_features_by_compound)
#' # @export
get_features_by_compound <- function(compound, organoid_type, 
                                     feature_type, configdir) {
  warning("FUNCTION UNIMPLEMENTED")
  return(NULL)
  source(file.path(configdir, "watchdogConfig.R"))
  
  organoid_pattern = NULL
  if(organoid_type == "human") {
    organoid_pattern = "D"
  } else if(organoid_type == "mouse") {
    organoid_pattern = "M"
  } else if(organoid_type == "both") {
    organoid_pattern = "[D,M]"
  } else {
    warning("Unrecognized organoid_type parameter")
    return(NULL)
  }
  
  all_plates = list.files(
    featuresdir, pattern = sprintf("^%s", organoid_pattern))
  
  # Go through all plates and extract the wells corresponding to the 
  # desired compound
  wells_to_retrieve = list()
  for(plate in all_plates) {
    layout_regex = regexec(pattern = "L[0-9]*", text = plate)[[1]]
    layout_id = substr(plate, layout_regex[1], 
                       layout_regex[1] + attr(layout_regex, "match.length"))
    # layout = read.table(file.path(layoutdir, paste0(layout_id, ".txt")), 
    #                     header = TRUE, quote = "", stringsAsFactors = FALSE, 
    #                     sep = ";", comment.char = "")
    layout = read.csv2(file.path(layoutdir, paste0(layout_id, ".txt")), 
                       stringsAsFactors = FALSE)
    rel_layout = layout[layout$Product.Name == compound,]
    if(nrow(rel_layout) == 0) next
    
    wells_to_retrieve[[plate]] = cbind(rel_layout, "Plate" = plate, 
                                       stringsAsFactors = FALSE)
  }
  
  if(length(wells_to_retrieve) == 0) return(data.frame())
  
  wells_to_retrieve = do.call(rbind, wells_to_retrieve)
  rownames(wells_to_retrieve) = NULL
  
  # Retrieve features for each entry
  features = list()
  for(n in seq_len(nrow(wells_to_retrieve))) {
    p = wells_to_retrieve[n,"Plate"]
    r = wells_to_retrieve[n,"Row_ID_384"]
    c = wells_to_retrieve[n,"Col_ID_384"]
    if(nchar(c) == 1) c = paste0("0", c)
    f = get_features_by_well(platename = p, row = r, col = c, 
                             feature_type = feature_type, 
                             configdir = configdir)
    features[[paste0(c(p, r, c), collapse="_")]] = f
  }
  
  features = do.call(rbind, features)
  return(features)
}

#' @title Get a list of all compounds
#' 
#' @description Get a list of all compounds
#' 
#' @param configdir The directory of the configuration files
#' 
#' @return A vector of strings
#' 
#' @author Jan Sauer
#' 
#' #@examples print(get_all_compounds)
#' #@export
get_all_compounds <- function(configdir) {
  warning("FUNCTION UNIMPLEMENTED")
  return(NULL)
  source(file.path(configdir, "watchdogConfig.R"))
  
  all_layouts = list.files(layoutdir, pattern = "L[0-9]*.txt")
  
  all_drugs = c()
  for(lf in all_layouts) {
    layout = read.table(file.path(layoutdir, lf), header = TRUE,
                        quote = "", stringsAsFactors = FALSE, sep = ";",
                        comment.char = "")
    all_drugs = c(all_drugs, unique(layout$Product.Name))
  }
  
  return(unique(all_drugs))
}

# ---
# DEV NOTE: This was used to get inception features. Keep this when reimplementing
# ---
# get_features_by_well <- function(platename, row, col, feature_type, configdir) {
#   source(file.path(configdir, "watchdogConfig.R"))
#   
#   if(nchar(col) == 1) col = paste0("0", col)
# 
#   if(feature_type == "inception") {
#     wells = h5read(file.path(featuresdir, platename, feature_type, sprintf(
#       "%s_features_inception.h5", platename)), "well_names")
#     well_index = which(wells == sprintf("%s_%s", row, col))
#     features = h5read(file.path(featuresdir, platename, feature_type, sprintf(
#       "%s_features_inception.h5", platename)), "features", 
#       index = list(NULL, well_index))
#     H5close()
#     return(as.vector(features))
#   } else if(feature_type == "intensity_distribution") {
#     wells = h5read(file.path(featuresdir, platename, feature_type, sprintf(
#       "%s_intensity_distribution.h5", platename)), "well_ids")
#     well_index = which(wells == sprintf("%s_%s_%s", platename, row, col))
#     features = h5read(file.path(featuresdir, platename, feature_type, sprintf(
#       "%s_intensity_distribution.h5", platename)), "features", 
#       index = list(NULL, NULL, NULL, well_index))
#     return(features[,,,1])
#     H5close()
#   }
# }