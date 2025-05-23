% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filenames.R
\name{segmentationHdf5Filename}
\alias{segmentationHdf5Filename}
\title{Create hdf5 filename from row and col for the DNN segmentation}
\usage{
segmentationHdf5Filename(filedir, platename, row, col)
}
\arguments{
\item{filedir}{The directory of the file}

\item{platename}{The plate that the well is on}

\item{row}{The row of the well}

\item{col}{The column of the well}
}
\value{
The HDF5 filename
}
\description{
Creates the hdf5 filename with given parameters for the DNN segmentation
}
\examples{
print(segmentationHdf5Filename)
}
\author{
Jan Sauer
}
