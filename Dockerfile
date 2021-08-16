FROM rocker/rstudio:4.0.0

# R packages
## CRAN
RUN R -e "install.packages('janitor',dependencies=TRUE, repos='http://cran.rstudio.com/')" 

## Bioconductor
RUN R -e "BiocManager::install('org.Hs.eg.db',dependencies=TRUE)"
RUN R -e "BiocManager::install('biomaRt',dependencies=TRUE)"
