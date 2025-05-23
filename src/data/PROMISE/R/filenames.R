
#' @title Get filenames associated with a specific well
#' 
#' @description This function finds all (existing) files associated with a specific well based on the number of z-stacks, channels, and fields defined in a configuration file
#' 
#' @param plateIndir The directory name of the plate as generated by the microscope
#' @param platename The ID (and folder name) of the plate
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory with configuration files
#' 
#' @return A string containing the list of existing files for the given well
#' 
#' @author Jan Sauer
#' 
#' @examples print(imageFilenames)
#' @export
imageFilenames <- function(plateIndir, platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  identifiers = expand.grid(fields, channels, stacks, stringsAsFactors = FALSE)
  all_filenames = apply(identifiers, 1, function(x) singleImageFilename(configdir, row, col, x[[1]], x[[2]], x[[3]]))
  all_filenames = file.path(indir, plateIndir, all_filenames)
#  files_exist = sapply(all_filenames, function(x) file.exists(file.path(indir, plateIndir, x)))
  
#  return(all_filenames[files_exist])
  return(all_filenames)
}

#' @title Create tif image filename name
#' 
#' @description Creates the tif image filename given with the parameters
#' 
#' @param configdir The directory containing watchdogConfig.R
#' @param row The row of the well
#' @param col The column of the well
#' @param field The field of the well
#' @param channel The color channel
#' @param stack The z-stack of the image
#' 
#' @return The image filename
#' 
#' @author Jan Sauer
#' 
#' @examples print(singleImageFilename)
#' @export
singleImageFilename <- function(configdir, row, col, field, channel, stack) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  # Match the excitement wavelength to the corresponding imaging channel
  exciteWavelength = exciteWavelengths[match(channel, channels)]
  filename_template = "%s - %s(fld %s wv %s - %s z %s).tif"
  return(sprintf(filename_template, row, col, field, exciteWavelength, channel, stack))
}

#' @title Create hdf5 filename from row and col
#' 
#' @description Creates the hdf5 filename with given parameters
#' 
#' @param filedir The directory of the file
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir This parameter is useless but still here for backwards compatibility
#' 
#' @return The HDF5 filename
#' 
#' @author Jan Sauer
#' 
#' @examples print(hdf5Filename)
#' @export
hdf5Filename <- function(filedir, platename, row, col, configdir = "NONE") {
  filename_template = "%s_%s_%s.h5"
  return(file.path(filedir, sprintf(filename_template, platename, row, col)))
}

#' @title Create hdf5 filename from row and col for the projections
#' 
#' @description Creates the hdf5 filename with given parameters for the contrast projections
#' 
#' @param filedir The directory of the file
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir This parameter is useless but still here for backwards compatibility
#' 
#' @return The HDF5 filename
#' 
#' @author Jan Sauer
#' 
#' @examples print(projectionHdf5Filename)
#' @export
projectionHdf5Filename <- function(filedir, platename, row, col, configdir = "NONE") {
    filename_template = "%s_%s_%s_contrastProjections.h5"
    return(file.path(filedir, sprintf(filename_template, platename, row, col)))
}

#' @title Create hdf5 filename from row and col for the validation
#' 
#' @description Creates the hdf5 filename with given parameters for the hdf5 validation
#' 
#' @param filedir The directory of the file
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir This parameter is useless but still here for backwards compatibility
#' 
#' @return The HDF5 filename
#' 
#' @author Jan Sauer
#' 
#' @examples print(validationHdf5Filename)
#' @export
validationHdf5Filename <- function(filedir, platename, row, col, configdir="NONE") {
  filename_template = "%s_%s_%s_isvalid.txt"
  return(file.path(filedir, sprintf(filename_template, platename, row, col)))
}

#' @title Create hdf5 filename from row and col for the intensity distribution
#' 
#' @description Creates the hdf5 filename with given parameters for the intensity distribution
#' 
#' @param filedir The directory of the file
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' 
#' @return The HDF5 filename
#' 
#' @author Jan Sauer
#' 
#' @examples print(intensityHdf5Filename)
#' @export
intensityHdf5Filename <- function(filedir, platename, row, col) {
  filename_template = "%s_%s_%s_intensity_distribution.h5"
  return(file.path(filedir, sprintf(filename_template, platename, row, col)))
}

#' @title Create image filenames from row and col for thumbnail images
#' 
#' @description Creates image filenames with given parameters for thumbnail images
#' 
#' @param filedir The directory of the file
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files. Unnecessary but here for backwards compatibility
#' @param level 1: The larger thumbnail; 2: The smaller thumbnail
#' @param addPath Logical, if TRUE, the full path is added to the filename
#' 
#' @return The image filename
#' 
#' @author Bernd Fischer
#' 
#' @examples print(thumbnailFilename)
#' @export
thumbnailFilename <- function(filedir, platename, row, col, configdir=NULL, level=1, addPath=TRUE) {
    filename_template = "%s_%s_%s_%d.jpeg"
    if (addPath) {
        res = file.path(filedir, sprintf(filename_template, platename, row, col, level))
    } else {
        res = sprintf(filename_template, platename, row, col, level)
    }
    res
}

#' @title Create hdf5 filename for features
#' 
#' @description Creates hdf5 filenames with given parameters for the organoid features
#' 
#' @param filedir The directory of the file
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#' 
#' @return The hdf5 filename
#' 
#' @author Jan Sauer
#' 
#' @examples print(featuresHdf5Filename)
#' @export
featuresHdf5Filename <- function(filedir, platename, row, col, configdir=NULL) {
  filename_template = "%s_%s_%s_features.h5"
  return(file.path(filedir, sprintf(filename_template, platename, row, col)))
}

#' @title Create hdf5 filename from row and col for the DNN segmentation
#' 
#' @description Creates the hdf5 filename with given parameters for the DNN segmentation
#' 
#' @param filedir The directory of the file
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' 
#' @return The HDF5 filename
#' 
#' @author Jan Sauer
#' 
#' @examples print(segmentationHdf5Filename)
#' @export
segmentationHdf5Filename <- function(filedir, platename, row, col) {
  filename_template = "%s_%s_%s_DNNsegmentation.h5"
  return(file.path(filedir, sprintf(filename_template, platename, row, col)))
}