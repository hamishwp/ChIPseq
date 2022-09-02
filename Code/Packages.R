###########################################################################
######################## INSTALL REQUIRED PACKAGES ########################
###########################################################################
if(installer){
  # # Now let's install all the other necessary packages:
  list.of.packages <- c("ggplot2","tidyverse","caret","stringr","reticulate",
                        "tensorflow","keras",
                        "FactoMineR","factoextra","gridExtra")
  # # Install the packages that aren't already installed
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  install_tensorflow(envname = "r-reticulate",version = "cpu")
  install_keras(envname = "r-reticulate",version = "cpu")
  virtualenv_create("r-reticulate", python = "/home/hamishwp/miniconda3/bin/python3.9")
  library(tensorflow)
  # install_tensorflow(envname = "r-reticulate")
  library(keras)
  # install_keras(envname = "r-reticulate")
  library(tensorflow)
  tf$constant("Hello Tensorflow!")
  
}
###########################################################################
############################# LOAD LIBRARIES ##############################
###########################################################################
library(stringr)
library(tidyverse)
library(magrittr)
library(reticulate)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(tensorflow)
library(keras)
library(caret)

