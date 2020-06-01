

 model.loop<-function(i) { 
  library(pscl)
  library(caret)
  library(randomForest)
  source("~/Documents/R/bin/OE.load.and.source.R")
  source("~/Documents/R/bin/OE.caret.load.and.source.R")

  ### run CARET to find optimal vars 
  
  foo<-which(is.na(obs.rc[,i]))
  if (length(foo) > 0) { sub1<-obs.rc[-foo,]}
  if (length(foo) == 0) { sub1<-obs.rc}
  
  if (length(sub1[,1]) > 0 & sum(sub1[,i])>0 ) { 
  
  sub2<-subset(stations.rc, row.names(stations.rc) %in% row.names(sub1))
  
  set.seed(10)
  rfe.out<-rfe(y=sub1[,i], x=sub2[,cand.var], 
             sizes=2:10,  rfeControl=ctrl, maximize=T, na.action=na.omit) 

# record optvars
optvars.rfe <- rfe.out$optVariables
optvars.rfe

### build RF with optimal vars from caret 
set.seed(10)
rf.out<-randomForest(y=sub1[,i], x=sub2[,optvars.rfe], 
                     ntree=500,importance=TRUE, norm.votes=TRUE, 
                     keep.forest=TRUE,na.action=na.omit, proximity=T)
saveas<-paste0(assemblage, ".rf.models/",assemblage, ".", i, "RF.model.Rdata")
save(rf.out, file=saveas)

rsq<-as.data.frame(rf.out$rsq)
rsq<-(rsq[500,]) 

# ### apply new RF model to CAL data 
# set.seed(10)
# predicted<-predict(rf.out, stations.rc[,optvars.rfe])
# observed<-obs.rc[,i]
# 
# set.seed(10)
# fit1<-glm(predicted~observed)
# pr2<- as.data.frame(pR2(fit1))
# set.seed(10)
# fit2<-lm(predicted~observed)
# r2 <- format(summary(fit2)$adj.r.squared, digits=3)
# plot(predicted,observed)

  #capture.output(paste(i, "rsq",rsq, "r2",r2, "pr2",pr2[4,], "pr2",pr2[6,]), file="refcal.lm.out.txt", append=T) 
  capture.output(paste(i, "rsq",rsq), file="refcal.lm.out.txt", append=T) 
  
}}
##### END LOOP #######
######################### convert txt to csv 