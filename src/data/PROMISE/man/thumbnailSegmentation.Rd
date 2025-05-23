% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thumbnails.R
\name{thumbnailSegmentation}
\alias{thumbnailSegmentation}
\title{Thumbnail images of DNN Segmentation Masks}
\usage{
thumbnailSegmentation(plateIndir, platename, row, col, configdir)
}
\arguments{
\item{plateIndir}{The base directory on the server, which contains a folder for each plate}

\item{platename}{The ID (and folder name) of the plate}

\item{row}{The row of the well}

\item{col}{The column of the well}

\item{configdir}{The configuration file, which defines the number of z-stacks, channels, and fields}
}
\value{
A boolean value indicating if the projection was successful
}
\description{
This function generates thumbnail images of the segmentation masks
}
\examples{
print(thumbnailSegmentation)
}
\author{
Jan Sauer
}
