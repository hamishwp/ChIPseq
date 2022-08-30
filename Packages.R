###########################################################################
######################## INSTALL REQUIRED PACKAGES ########################
###########################################################################
# List of all the packages that we require, starting with BiocManager for the ChIPseeker package:
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ChIPseeker")
# Now let's install all the other necessary packages:
list.of.packages <- c("ggplot2","tidyverse","TxDb.Hsapiens.UCSC.hg19.knownGene",
                      "clusterProfiler") 
# Install the packages that aren't already installed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)



if(length(list.of.packages[!("ggmap" %in% installed.packages()[,"Package"])])){devtools::install_github("dkahle/ggmap")}


###########################################################################
############################# LOAD LIBRARIES ##############################
###########################################################################

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)


