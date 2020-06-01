# wrapper script for creating MMIs 

setwd("~/Documents/R/ASCI/pMMI_newnew/")

diatoms = F # one at a time
sba = F
hybrid = T

{ 
###Loading Needed Packages
library(vegan)
library(indicspecies)
library(reshape2)
library(tidyverse)
library(plyr) # WARNING: problems arise if you have dplyr already loaded, so maybe start from scratch
options(gsubfn.engine = "R")
library(sqldf)
library(ggplot2)
library(pscl)
library(caret)
library(randomForest)
library(ggplot2)
source("~/Documents/R/bin/OE.load.and.source.R")
source("~/Documents/R/bin/OE.caret.load.and.source.R")
source("~/Documents/R/bin/cbind.fill.R")
source("~/Documents/R/bin/cbind.na.R")


if (diatoms == T) { assemblage <- "diatoms" }
if (sba == T) { assemblage <- "sba" }
if (hybrid == T) { assemblage <- "hybrid" }

if (diatoms==T) {dir.create("~/Documents/R/ASCI/pMMI_newnew/diatoms")}
if (sba==T) {dir.create("~/Documents/R/ASCI/pMMI_newnew/sba")}
if (hybrid==T) {dir.create("~/Documents/R/ASCI/pMMI_newnew/hybrid")}

if (diatoms==T) {dir.create("~/Documents/R/ASCI/pMMI_newnew/diatoms/diatoms.plots")}
if (sba==T) {dir.create("~/Documents/R/ASCI/pMMI_newnew/sba/sba.plots")}
if (hybrid==T) {dir.create("~/Documents/R/ASCI/pMMI_newnew/hybrid/hybrid.plots")}
}


# Step 1 
source("~/Documents/R/ASCI/pMMI_newnew/Step1_EvalMetrics_4.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 2a 
source("~/Documents/R/ASCI/pMMI_newnew/Step2a_ModelMetrics3.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 2b 
source("~/Documents/R/ASCI/pMMI_newnew/Step2b_ModelMetrics3.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 3
source("~/Documents/R/ASCI/pMMI_newnew/Step3_ScoreMetrics.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 4a 
source("~/Documents/R/ASCI/pMMI_newnew/Step4a_ScreenMetrics_6.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 4b
source("~/Documents/R/ASCI/pMMI_newnew/Step4b_Compare_pMMIs_4.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 4bb # requires hand-holding 
source("Step4bb_iterate.metrics_all.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 4bb2 save winning metrics
source("Step4bb2.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 4c
source("~/Documents/R/ASCI/pMMI_newnew/Step4c_Winning_PMMIs.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# following steps are done on all assemblages at the same time 

# Step 5 
source("Step5_ASCI_Scores_3.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 6
source("Step6_PullRFmodels.R")
setwd("~/Documents/R/ASCI/pMMI_newnew/")

################################
# 

# Performance 1 
source("~/Documents/R/ASCI/PERFORMANCE_newnew/001_Performance14.R")
source("~/Documents/R/ASCI/PERFORMANCE_newnew/002_Figures10.R")
source("~/Documents/R/ASCI/PERFORMANCE_newnew/003_CompareCSCI.R")
source("~/Documents/R/ASCI/PERFORMANCE_newnew/004_Fig5_LowE2.R")



