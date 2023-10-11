# Start from an R base image
FROM rocker/r-ver:4.2.0

# Install system dependencies for R packages and devtools
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libz-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev 


# Install the devtools package with verbose output
RUN R -e "options(warn = 2); install.packages('devtools', repos='https://cloud.r-project.org/', verbose=TRUE)"

# Install R libraries from DESCRIPTION file
RUN R -e "install.packages(c('Biostrings', 'data.table', 'dplyr', 'GenomicRanges', 'IRanges', 'Matrix', 'matrixStats', 'pbmcapply', 'Rsamtools', 'Rsubread', 'S4Vectors', 'utils', 'BiocManager'), repos='https://cloud.r-project.org/')"
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("Biostrings", "GenomicRanges", "Rsamtools", "Rsubread"))'
# Install additional R libraries from README.md
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("WES.1KG.WUGSC", "BSgenome.Hsapiens.UCSC.hg19"))'#
RUN apt-get install -y libbz2-dev
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c( "BSgenome.Hsapiens.UCSC.hg19"))'

# Clone and install HMZDupFinder and SLMSeg from GitHub
RUN R -e "devtools::install_github('cluhaowie/HMZDupFinder')"#
RUN R -e "devtools::install_github('cluhaowie/VizCNV/SLMSeg')"
