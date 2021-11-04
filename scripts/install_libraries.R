### INSTALL LIBRARY PACKAGES AT THE BEGINNING OF R (only first time)

# install bioMart for bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# install tidyverse
install.packages("tidyverse")

# install valr for bed_merge()
install.packages("valr")

# update R
install.packages("installr", dependencies = TRUE)
library(installr)
updateR()

# install xlsx to save in excel
install.packages("xlsx")
