#find best custom combination of metrics for pMMI 
#step 1:  import your data from Step4b so it is in your global environment 

setwd(paste0("~/Documents/R/ASCI/pMMI_newnew/", assemblage, "/"))

### Winning metrics --------------------------------------------------------------

if (assemblage == "diatoms") { 


  win.metrics<-c("Salinity.BF.richness",
                 "cnt.spp.most.tol",
                 "prop.spp.Trophic.E",
                 "prop.spp.IndicatorClass_TN_low",
                 "prop.spp.Planktonic",
                 "EpiRho.richness",
               NULL ) 
  write.csv(win.metrics, "win.metrics.csv")

  win.metrics.raw<-c("Salinity.BF.richness_raw",
                     "cnt.spp.most.tol_raw",
                     "prop.spp.Trophic.E_raw",
                     "prop.spp.IndicatorClass_TN_low_raw",
                     "prop.spp.Planktonic_raw",
                     "EpiRho.richness_raw",
                     NULL ) 
  write.csv(win.metrics.raw, "win.metrics.raw.csv") }

if (assemblage == "sba") { 
  
  win.metrics<-c("prop.spp.IndicatorClass_TP_high_raw",
                 "prop.spp.IndicatorClass_DOC_high_raw",
                 "prop.spp.IndicatorClass_NonRef_raw",
                 "prop.spp.ZHR_raw",
                 NULL ) 
  write.csv(win.metrics, "win.metrics.csv")
  
  win.metrics.raw<-c("prop.spp.IndicatorClass_TP_high_raw",
                     "prop.spp.IndicatorClass_DOC_high_raw",
                     "prop.spp.IndicatorClass_NonRef_raw",
                     "prop.spp.ZHR_raw",
                     NULL ) 
  write.csv(win.metrics.raw, "win.metrics.raw.csv") }

if (assemblage == "hybrid") { 
  
  win.metrics<-c("cnt.spp.most.tol",
                  "prop.spp.ZHR_raw",
                 "Salinity.BF.richness",
                 "prop.spp.Planktonic",
                 "cnt.spp.IndicatorClass_TP_high",
                 "prop.spp.Trophic.E",
                 "EpiRho.richness",
                 "OxyRed.DO_30.richness",
                 NULL ) 
  write.csv(win.metrics, "win.metrics.csv")
  
  win.metrics.raw<-c(
                    "cnt.spp.most.tol_raw",
                      "prop.spp.ZHR_raw",
                     "Salinity.BF.richness_raw",
                     "prop.spp.Planktonic_raw",
                     "cnt.spp.IndicatorClass_TP_high_raw",
                     "prop.spp.Trophic.E_raw",
                     "EpiRho.richness_raw",
                     "OxyRed.DO_30.richness_raw",
                     NULL ) 
  write.csv(win.metrics.raw, "win.metrics.raw.csv") }


