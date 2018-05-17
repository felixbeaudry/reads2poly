library(data.table)
library(devtools)
library(ggfortify)
library(ape)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(sqldf)

setwd('~/Google Drive/Research/Data/')

####################GBS
localdf <- fread("GBSlocations.txt")

GBS <- fread("GBS.spotless.012_popd.txt")
GBSID<-fread('gbsnames.txt',header=T)
#GBSs <- data.frame(c(GBSID,GBS[,-1]), row.names=1)
GBSs <- data.frame(c(GBSIDXYY[-c(53)],GBS[-c(53),-1]), row.names=1) #OKRAT_17 
#GBSs <- data.frame(c(GBSID[-c(30:32,46:51,53)],GBS[-c(30:32,46:51,53),-1]), row.names=1) #OKRAT_17 + LABEN + OKBAC

GBSIDXYY<-fread('gbsnames2.txt',header=T)
GBSXYY <- data.frame(c(GBSIDXYY[-c(30:32,46:57,76:92)],GBS[-c(30:32,46:57,76:92),-1]), row.names=1) #XY removed

pwGBS<-dist.gene(GBS,method = "pairwise")
GBStree<-nj(as.dist(pwGBS))
GBStree$tip.label <- GBSID$Name
plot(GBStree,show.tip.label=TRUE,type = "unrooted", lab4ut="axial",cex=0.9,rotate.tree=100)
write.tree(GBStree, file = "GBStree")
         #  , append = FALSE, digits = 10, tree.names = FALSE)

GBSpca <- PCA(GBSs[,-(1:4)], graph = FALSE)
#fviz_pca_ind(GBSpca)
#fviz_pca_ind(GBSpca, label="none")
#fviz_pca_ind(GBSpca, label="GBS Population PCA",habillage=as.factor(GBSs$Pop))
#fviz_pca_ind(GBSpca, label="", habillage=as.factor(GBSs$State),title="",addEllipses=TRUE, ellipse.level=0.95, ggtheme =  theme_classic())

 GBScoord <- GBSpca$ind$coord
 GBScoords <- setDT(data.frame(GBScoord), keep.rownames = TRUE)[]
 GBScoordsN <- separate(GBScoords, rn, c("pop","ind"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
 GBScoordsN$state <- GBSs$State
 
 GBSeig <- GBSpca$eig
 
XYXYYPCAplot <- 
   ggplot(GBScoordsN,aes(x=-(Dim.1), y=Dim.2, label=pop,color=state)) + geom_text(size=10) + 
   theme_bw(base_size = 18) + guides(color = FALSE) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "PCA1 (4.98%)",y = "PCA2 (2.74%)")
 

####XYY####
 GBSXYYpca <- PCA(GBSXYY[,-(1:4)], graph = FALSE)
 #fviz_pca_ind(GBSXYYpca, label="", habillage=as.factor(GBSXYY$Sex),title="",addEllipses=TRUE, ellipse.level=0.95, ggtheme =  theme_classic())
 #fviz_pca_ind(GBSXYYpca, label="", habillage=as.factor(GBSXYY$State), title="",ellipse.level=0.5, addEllipses=TRUE, ggtheme =  theme_classic())
 
 GBSXYYcoord <- GBSXYYpca$ind$coord
 GBSXYYcoords <- setDT(data.frame(GBSXYYcoord), keep.rownames = TRUE)[]
 GBSXYYcoordsN <- separate(GBSXYYcoords, rn, c("pop","ind"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
 GBSXYYcoordsN$state <- GBSXYY$State
 
 GBSXYYeig <- GBSXYYpca$eig
 
 #XYYPCAplot <- 
   ggplot(GBSXYYcoordsN,aes(x=Dim.2, y=Dim.1, label=pop,color=state)) + geom_text(size=10) + 
 theme_bw(base_size = 18) + guides(color = FALSE) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "PCA2 (3.17%)",y = "PCA1 (3.31%)")

   locDimXYY <- data.frame(sqldf('select localdf.Latitude, localdf.Longitude, GBSXYYcoordsN.* 
                           from GBSXYYcoordsN 
                              left join localdf on GBSXYYcoordsN.pop = localdf.name'))
   
   ggplot(locDimXYY,aes(x=Latitude,y=Dim.1)) + geom_point() + geom_smooth(method="lm",se = TRUE,size=1)
   XYY_lat_1_model <- lm(Latitude ~ Dim.1 , data=locDimXYY)
   summary(XYY_lat_1_model)
   
   ggplot(locDimXYY,aes(x=Longitude,y=Dim.2)) + geom_point() + geom_smooth(method="lm",se = TRUE,size=1)
   XYY_long_2_model <- lm(Longitude ~ Dim.2 , data=locDimXYY)
   summary(XYY_long_2_model)
   
####XY####
 GBSXY <- data.frame(c(GBSIDXYY[c(30:32,46:52,54:57,76:92)],GBS[c(30:32,46:52,54:57,76:92),-1]), row.names=1) #XY removed
 GBSXYpca <- PCA(GBSXY[,-(1:4)], graph = FALSE)
 fviz_pca_ind(GBSXYpca, label="", habillage=as.factor(GBSXY$State), title="", addEllipses=TRUE, ggtheme =  theme_classic())
 
 GBSXYcoord <- GBSXYpca$ind$coord
 GBSXYcoords <- setDT(data.frame(GBSXYcoord), keep.rownames = TRUE)[]
 GBSXYcoordsN <- separate(GBSXYcoords, rn, c("pop","ind"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
 GBSXYcoordsN$state <- GBSXY$State
 
 GBSXYeig <- GBSXYpca$eig

 
 #XYPCAplot <- 
   ggplot(GBSXYcoordsN,aes(x=Dim.1, y=Dim.2, label=pop,color=state)) + geom_text(size=10) + 
   theme_bw(base_size = 18) + guides(color = FALSE) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "PCA1 (5.95%)",y = "PCA2 (5.03%)")
 
 locDimXY <- data.frame(sqldf('select localdf.Latitude, localdf.Longitude, GBSXYcoordsN.* 
                           from GBSXYcoordsN 
                               left join localdf on GBSXYcoordsN.pop = localdf.name'))
 
 ggplot(locDimXY,aes(x=Latitude,y=Dim.2)) + geom_point() + geom_smooth(method="lm",se = TRUE,size=1)
 XY_lat_2_model <- lm(Latitude ~ Dim.2 , data=locDimXY)
 summary(XY_lat_2_model)
 
 ggplot(locDimXY,aes(x=Longitude,y=Dim.1)) + geom_point() + geom_smooth(method="lm",se = TRUE,size=1)
 XY_long_1_model <- lm(Longitude ~ Dim.1 , data=locDimXY)
 summary(XY_long_1_model)
 
 
 
 multiplot(XYXYYPCAplot,XYPCAplot,XYYPCAplot,cols=3)
 
 
#,addEllipses=TRUE, ellipse.level=0.95

GBSscale <- scale(GBSs[,-(1:3)])
GBSscale[is.na(GBSscale)] <- 1 
fviz_nbclust(GBSscale, kmeans, method = "gap_stat")

GBScut <- hcut(GBSs[,-(1:3)], k = 2, stand = TRUE)
fviz_dend(GBScut, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))


####################RNA

RNA <-  fread("allpops_aa.vclean.n012")
RNAID<-fread('name_aa.txt',header=T)
RNAs <- data.frame(c(RNAID,RNA[,-1]), row.names=1)
#RNAs <- data.frame(c(RNA[,-1]))

pwRNA<-dist.gene(RNAs,method = "pairwise")
RNAtree<-nj(as.dist(pwRNA))
RNAtree$tip.label <- RNAID$ID
plot(RNAtree,show.tip.label=TRUE,type = "unrooted", lab4ut="axial",cex=1.5,rotate.tree=-20)

RNApca <- PCA(RNAs[,-(1:3)], graph = FALSE)
