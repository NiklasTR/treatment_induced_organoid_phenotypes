FROM rocker/rstudio:4.0.0

RUN apt-get update && apt-get install -y \
  build-essential \
  libglpk40

# R packages
## CRAN
RUN R -e "install.packages('janitor',dependencies=TRUE)" 
RUN R -e "install.packages(c('igraph', 'ggdag'), dependencies = TRUE)"

## Bioconductor
RUN R -e "BiocManager::install('org.Hs.eg.db',dependencies=TRUE)"
RUN R -e "BiocManager::install('biomaRt',dependencies=TRUE)"
RUN R -e "BiocManager::install('enrichplot',dependencies=TRUE)"
RUN R -e "BiocManager::install('clusterProfiler',dependencies=TRUE)"
RUN R -e "BiocManager::install('ReactomePA',dependencies=TRUE)"
