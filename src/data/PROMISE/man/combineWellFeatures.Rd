% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featureExtraction.R
\name{combineWellFeatures}
\alias{combineWellFeatures}
\title{Combine well features}
\usage{
combineWellFeatures(platename, configdir)
}
\arguments{
\item{platename}{The plate that the well is on}

\item{configdir}{The directory of the configuration files}
}
\value{
A boolean value indicating if the features were correctly combined 
and written to disk
}
\description{
Combines the features for each individual well of a plate. 
This will only combine features if the features for every well have been 
calculated. The output is written to a file.
}
\examples{
print(combineWellFeatures)
}
\author{
Jan Sauer
}
