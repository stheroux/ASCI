# modeling metrics for reference sites 

setwd(paste0("~/Documents/R/ASCI/pMMI_newnew/", assemblage, "/"))

results<-read.table("refcal.lm.out.txt")
foo<-as.matrix(results)
foo2<-do.call(rbind, strsplit(as.character(foo[,2]), " "))
foo2<-as.data.frame(foo2[,c(1,3)])
names(foo2)<-c('top.metrics', "rsq")

foo2$rsq<-as.numeric(as.character(foo2$rsq))

if (diatoms ==T) { 
write.csv(foo2, 'diatoms.refcal.rsq.csv') }

if (sba ==T) { 
  write.csv(foo2, 'sba.refcal.rsq.csv') }

if (hybrid ==T) { 
  write.csv(foo2, 'hybrid.refcal.rsq.csv') }

metrics<-foo2$top.metrics
metrics<-droplevels(metrics)
#results$rsq<-as.numeric(as.character(results$rsq))
tops <-droplevels(subset(foo2, rsq>0.2)) # or pr2.1>0.1
top.metrics<-tops$top.metrics


if (diatoms == T ) {write.csv(top.metrics, "diatoms.top.metrics.csv") }
if (sba == T ) {write.csv(top.metrics, "sba.top.metrics.csv") }
if (hybrid == T ) {write.csv(top.metrics, "hybrid.top.metrics.csv") }

##########################################################
#using RF models, predict values 
##########################################################

predicted<-data.frame(row.names(stations))
row.names(predicted)<-predicted[,1]

for (i in unique(metrics)) {        # i<-"prop.spp.BCG45"
 
  # import RF model 
  filename<-paste0("~/Documents/R/ASCI/pMMI_newnew/", assemblage,"/",assemblage,".rf.models/",assemblage,".",i,"RF.model.Rdata" )
  load(filename)
  pred.var<-row.names(rf.out$importance)
  
  ### apply new RF model to ALL data 
  stations1<-stations[,pred.var]
  foo<-which(is.na(rowSums(stations1)))
  if (length(foo) > 0) {stations1<-stations1[-foo,]}
  set.seed(10)
  pred<-as.data.frame(predict(rf.out, stations1))
  names(pred)[1]<-i
  predicted<-cbind.na(predicted, pred) } 

predicted<-predicted[,-1]
colnames(predicted)<-metrics

if (diatoms == T ) {write.csv(predicted, "diatoms.predicted.metrics.csv") }
if (sba == T ) {write.csv(predicted, "sba.predicted.metrics.csv") }
if (hybrid == T ) {write.csv(predicted, "hybrid.predicted.metrics.csv") }

predicted.top<-predicted[,top.metrics]

if (diatoms == T ) {write.csv(predicted.top, "diatoms.predicted.top.metrics.csv") }
if (sba == T ) {write.csv(predicted.top, "sba.predicted.top.metrics.csv") }
if (hybrid == T ) {write.csv(predicted.top, "hybrid.predicted.top.metrics.csv") }






