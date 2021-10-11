FROM rocker/verse:4.0.0

RUN apt-get update && apt-get install -y \
  build-essential \
  libglpk40 \
  tmux \
  libxml2-dev \
  libcurl4-openssl-dev \
  libssl-dev \ 
  nano \
  less

# R packages
## CRAN
#RUN R -e "install.packages('tidyverse',dependencies=TRUE)" 
RUN R -e "install.packages('igraph',dependencies=TRUE)" 
### grouped packages where I do not expect errors
RUN R -e "install.packages(c('janitor', 'ggdag', 'here', 'cowplot', 'gridExtra', 'pheatmap', 'princurve', 'scico'), dependencies = TRUE)"
RUN R -e "install.packages(c('ggrastr'), dependencies = TRUE)"

## Bioconductor
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('org.Hs.eg.db',dependencies=TRUE)"
RUN R -e "BiocManager::install('biomaRt',dependencies=TRUE)"
RUN R -e "BiocManager::install('enrichplot',dependencies=TRUE)"
RUN R -e "BiocManager::install('clusterProfiler',dependencies=TRUE)"
RUN R -e "BiocManager::install('ReactomePA',dependencies=TRUE)"
### grouped packages where I do not expect errors
RUN R -e "BiocManager::install(c('lumi', 'limma', 'affy'),dependencies=TRUE)"

## devtools and dependencies
### installing dependencies for ggrastr
RUN apt-get install -y \
  libgtk2.0-dev \
  libcairo2-dev \
   xvfb \
   xauth \
   xfonts-base \
   libxt-dev
RUN R -e "devtools::install_github('VPetukhov/ggrastr', build_vignettes = FALSE)"
### installing dependencies for MOFA2
RUN R -e "devtools::install_version('matrixStats', version = '0.60.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('MOFA2',dependencies=TRUE)"
