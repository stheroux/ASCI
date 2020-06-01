# learning lapply and making a screen.metrics function

### trying in parallel 

#lapply(1:3, function(x) c(x, x^2, x^3))

screen.metrics<-function(i) {
  win.metrics<-i
  combined.win.metrics<- data.frame(combined.sc[,win.metrics])
  scores<-data.frame(stations$SiteSetSample2, rowMeans(combined.win.metrics))
  colnames(scores)<-c("Type", "Means")
  fit<-aov(scores$Means~scores$Type)
  Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
  Fv<-summary(fit)[[1]][["F value"]][[1]]
  saveas<-paste(colnames(combined.win.metrics), collapse = " ")
  Pv.list[saveas]<-Pv
  Fv.list[saveas]<-Fv
  Fv
  
  # sd refcal
  scores.rc<-subset(scores, row.names(scores) %in% row.names(rc))
  sd1<-sd(scores.rc$Means, na.rm = T)
  sd1.list[saveas]<-sd1
  sd1
  
  # PSA anova 
  stations.rc<-subset(stations, row.names(stations) %in% row.names(scores.rc))
  fit<-aov(scores.rc$Means~stations.rc$PSA6c)
  Fv2<-summary(fit)[[1]][["F value"]][[1]]
  psa.anova.list[saveas]<-Fv2
  Fv2

  capture.output(paste(toString(i),Fv,sd1,Fv2), file="out.txt", append=T)
}



#sapply(foo, function(x) fun=screen.metrics(x))
#lapply(foo, function(x) fun=screen.metrics(x))
