% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow.R
\name{PROMISEsegmentation}
\alias{PROMISEsegmentation}
\title{Segmentation Workflow}
\usage{
PROMISEsegmentation(plateIndir, platename, row, col, configdir)
}
\arguments{
\item{plateIndir}{The directory name of the plate as generated by the microscope}

\item{platename}{The ID (and folder name) of the plate}

\item{row}{The row of the well}

\item{col}{The column of the well}

\item{configdir}{The configuration file, which defines the number of z-stacks, channels, and fields}
}
\description{
This is the PROMISE workflow to segment HDF5 projections
}
\examples{
print(PROMISEsegmentation)
}
\author{
Jan Sauer
}
