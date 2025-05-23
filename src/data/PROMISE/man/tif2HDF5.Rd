% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tif2HDF5.R
\name{tif2HDF5}
\alias{tif2HDF5}
\title{Turn a series of TIF images into an HDF5 file}
\usage{
tif2HDF5(plateIndir, platename, row, col, configdir)
}
\arguments{
\item{plateIndir}{The directory name of the plate as generated by the microscope}

\item{platename}{The ID (and folder name) of the plate}

\item{row}{The row of the well}

\item{col}{The column of the well}

\item{configdir}{The configuration file, which defines the number of z-stacks, channels, and fields}
}
\value{
A boolean value indicating if the save was successful
}
\description{
This function stores all images associated with a well in an HDF5 file
}
\examples{
print(tif2HDF5)
}
\author{
Jan Sauer
}
