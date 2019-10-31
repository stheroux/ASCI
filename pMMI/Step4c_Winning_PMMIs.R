# Step 4c Compile winning metrics 

setwd("~/Documents/R/ASCI/pMMI/")

diatoms = F # pick one 
sba = F
hybrid = T

if (diatoms == T) {assemblage <- "diatoms"}
if (sba == T) {assemblage <- "sba"}
if (hybrid == T) {assemblage <- "hybrid"}

if (diatoms == T) {setwd("~/Documents/R/ASCI/pMMI/diatoms/")}
if (sba == T) {setwd("~/Documents/R/ASCI/pMMI/sba/")}
if (hybrid == T) {setwd("~/Documents/R/ASCI/pMMI/hybrid/")}

#################################
# Import data ------------------
#################################

importfile<-paste0(assemblage, ".combined.scaled.csv")
combined.sc<-read.csv(importfile, header=T, row.names=1)
metrics.list<-colnames(combined.sc)

# Import stations data 
stations<-read.csv('~/Documents/R/ASCI/DATA/algae.site.data.07202019.csv')
row.names(stations)<-stations$SampleID_old
stations<-stations[order(row.names(stations)),]
stations<-subset(stations, row.names(stations) %in% row.names(combined.sc), drop=T)
foo<-sapply(stations, is.numeric)
stations.num<-stations[,foo] # only numeric columns 

calval<-subset(stations, stations$RefCalVal==c("Cal", "Val"))
calval<-calval[order(row.names(calval)),]

row.names(stations) == row.names(combined.sc)

stations.ref<-subset(stations, row.names(stations) %in% row.names(calval), drop=T)
stations.str<-subset(stations, stations$SiteSetSample2=="Stressed")

combined.sc.ref<-subset(combined.sc, row.names(combined.sc) %in% row.names(calval), drop=T)
combined.sc.str<-subset(combined.sc, stations$SiteSetSample2=="Stressed")
  
############################ 
# winning pMMIs ###########
############################ 

  win.metrics<-read.csv("win.metrics.csv", row.names=1)
  win.metrics<-win.metrics$x
  
  win.metrics.null<-read.csv("win.metrics.raw.csv", row.names=1)
  win.metrics.null<-win.metrics.null$x


# Ref/Int/Str ANOVA 
foo<-which(colnames(combined.sc) %in% win.metrics)
combined.win.metrics<-as.data.frame(combined.sc[,foo])
scores<-data.frame(stations$SiteSetSample2, rowMeans(combined.win.metrics, na.rm = T))
colnames(scores)<-c("Type", "Means")
fit<-aov(scores$Means~scores$Type)
Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
Fv<-summary(fit)[[1]][["F value"]][[1]]
Fv

# PSA region ANOVA cal
cal<-subset(calval, calval$RefCalVal=="Cal")
scores.rc<-subset(scores, row.names(scores) %in% row.names(cal))
stations.rc<-subset(stations, row.names(stations) %in% row.names(scores.rc))
fit<-aov(scores.rc$Means~stations.rc$PSA6c)
Fv<-summary(fit)[[1]][["F value"]][[1]]
Fv

# val
val<-subset(calval, calval$RefCalVal=="Val")
scores.rv<-subset(scores, row.names(scores) %in% row.names(val))
stations.rv<-subset(stations, row.names(stations) %in% row.names(scores.rv))
fit<-aov(scores.rv$Means~stations.rv$PSA6c)
Fv<-summary(fit)[[1]][["F value"]][[1]]
Fv


if (diatoms==T) { write.csv(combined.win.metrics, "diatoms.combined.win.metrics.csv")}
if (sba==T) { write.csv(combined.win.metrics, "sba.combined.win.metrics.csv") }
if (hybrid==T) {  write.csv(combined.win.metrics, "hybrid.combined.win.metrics.csv")}

if (diatoms==T) { write.csv(scores, "diatoms.combined.win.metrics.scores.csv")}
if (sba==T) { write.csv(scores, "sba.combined.win.metrics.scores.csv") }
if (hybrid==T) {  write.csv(scores, "hybrid.combined.win.metrics.scores.csv")}

# Null metrics 
foo<-which(colnames(combined.sc) %in% win.metrics.null)
combined.win.metrics.null<-as.data.frame(combined.sc[,foo])
scores.null<-data.frame(stations$SiteSetSample2, rowMeans(combined.win.metrics.null, na.rm = T))
colnames(scores.null)<-c("Type", "Means")
fit<-aov(scores.null$Means~scores.null$Type)
Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
Fv<-summary(fit)[[1]][["F value"]][[1]]
Fv

if (diatoms==T) { write.csv(combined.win.metrics.null, "diatoms.combined.win.metrics.null.csv")}
if (sba==T) { write.csv(combined.win.metrics.null, "sba.combined.win.metrics.null.csv") }
if (hybrid==T) {  write.csv(combined.win.metrics.null, "hybrid.combined.win.metrics.null.csv")}

if (diatoms==T) { write.csv(scores.null, "diatoms.combined.win.metrics.scores.null.csv")}
if (sba==T) { write.csv(scores.null, "sba.combined.win.metrics.scores.null.csv") }
if (hybrid==T) {  write.csv(scores.null, "hybrid.combined.win.metrics.scores.null.csv")}


#############################################
# winning pMMI graphs 
#############################################

ggplot(scores, aes(y=as.numeric(scores$Means))) + 
  geom_boxplot(aes(x=factor(scores$Type), fill=factor(scores$Type))) + geom_jitter(aes(x=factor(scores$Type)))+
 scale_fill_manual(values=c("#FFCC66", "#339966", "#CC0000" )) +
  scale_x_discrete(limits=c("Reference", "Intermediate", "Stressed"))  +
  ylab("ScoreMeans") + xlab("")  + theme_bw(base_size = 15) + theme(legend.position="none")
saveas<-paste0(assemblage,".ScoreMeans.pdf")
ggsave(saveas, width = 8, height = 8)

######################################################
# compare pMMI score to env gradients 
######################################################

list3<- c( 
  "DayOfYear", "Elevation", "XSLOPE", "NHDSLOPE", "New_Lat", "New_Long", "LST32AVE", 
  "MEANP_WS", "PPT_00_09", "CondQR50", "TMAX_WS", "XWD_WS", 
  "KFCT_AVE", "BDH_AVE", "PRMH_AVE", "CondQR50", "LogWSA", "AREA_SQKM")

r2.list<-list()
slope.list<-list()

x<-scores$Means

for (y in list3) {  
  fit<-lm(stations[[y]]~x, na.action = "na.exclude")
  r2 <- format(summary(fit)$adj.r.squared, digits=3) 
  slope <-format(coef(fit)[2], digits=3)
  saveas<-paste(y)
  r2.list[[saveas]]<-r2
  slope.list[[saveas]]<-slope
}

lm.out<-cbind(r2.list, slope.list)
if (diatoms==T) {write.csv(lm.out, "diatoms.pMMI.lm.out.env.gradients.csv")}
if (sba==T) {write.csv(lm.out, "sba.pMMI.lm.out.env.gradients")}
if (hybrid==T) {write.csv(lm.out, "hybrid.pMMI.lm.out.env.gradients.csv")}



######################################################
# compare pMMI score to disturbance gradients 
######################################################

list3<- c( 
  "Nitrate_as_N_mgPerL","Total_Organic_Carbon_mgPerL" )

r2.list<-list()
slope.list<-list()

x<-scores$Means

for (y in list3) {  
  fit<-lm(stations[[y]]~x, na.action = "na.exclude")
  r2 <- format(summary(fit)$adj.r.squared, digits=3) 
  slope <-format(coef(fit)[2], digits=3)
  saveas<-paste(y)
  r2.list[[saveas]]<-r2
  slope.list[[saveas]]<-slope
}

lm.out<-cbind(r2.list, slope.list)
if (diatoms==T) {write.csv(lm.out, "diatoms.pMMI.lm.out.dist.gradients.csv")}
if (sba==T) {write.csv(lm.out, "sba.pMMI.lm.out.dist.gradients")}
if (hybrid==T) {write.csv(lm.out, "hybrid.pMMI.lm.out.dist.gradients.csv")}

ggplot(scores) + geom_point(aes(x=log10(stations$Nitrate_as_N_mgPerL), y=scores$Means)) + 
  geom_smooth(aes(x=log10(stations$Nitrate_as_N_mgPerL), y=scores$Means), method="lm")

ggplot(scores) + geom_point(aes(x=log10(stations$Total_Organic_Carbon_mgPerL), y=scores$Means)) + 
  geom_smooth(aes(x=log10(stations$Total_Organic_Carbon_mgPerL), y=scores$Means), method="lm")




# Metric correlations to each other ----------------------------

#install.packages("Hmisc")
# library(Hmisc)
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }
# 
# cor1<-cor(combined.sc[,win.metrics], method='pearson') #quick pearson correl
# 
# cor2 <- rcorr(as.matrix(combined.sc[,win.metrics])) #gives correlations and p vals
# cor2$r
# cor2$P
# 
# cor2.flat<-flattenCorrMatrix(cor2$r, cor2$P)
# saveas<-paste0(assemblage, ".core2.flat.csv")
# write.csv(cor2.flat, saveas)

