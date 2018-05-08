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
fviz_pca_ind(GBSpca)
fviz_pca_ind(GBSpca, label="none")
fviz_pca_ind(GBSpca, label="GBS Population PCA",habillage=as.factor(GBSs$Pop))
 fviz_pca_ind(GBSpca, label="", habillage=as.factor(GBSs$State),title="",addEllipses=TRUE, ellipse.level=0.95, ggtheme =  theme_classic())

 GBScoord <- GBSpca$ind$coord
 GBScoords <- setDT(data.frame(GBScoord), keep.rownames = TRUE)[]
 GBScoordsN <- separate(GBScoords, rn, c("pop","ind"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
 GBScoordsN$state <- GBSs$State
 
 GBSeig <- GBSpca$eig
 
 XYXYYPCAplot <- 
   ggplot(GBScoordsN,aes(x=-(Dim.1), y=Dim.2, label=pop,color=state)) + geom_text(size=10) + 
   theme_bw(base_size = 18) + guides(color = FALSE) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "PCA1 (4.98%)",y = "PCA2 (2.74%)")
 

 
 #XYY
 GBSXYYpca <- PCA(GBSXYY[,-(1:4)], graph = FALSE)
 #fviz_pca_ind(GBSXYYpca, label="", habillage=as.factor(GBSXYY$Sex),title="",addEllipses=TRUE, ellipse.level=0.95, ggtheme =  theme_classic())
 fviz_pca_ind(GBSXYYpca, label="", habillage=as.factor(GBSXYY$State), title="",ellipse.level=0.5, addEllipses=TRUE, ggtheme =  theme_classic())
 
 GBSXYYcoord <- GBSXYYpca$ind$coord
 GBSXYYcoords <- setDT(data.frame(GBSXYYcoord), keep.rownames = TRUE)[]
 GBSXYYcoordsN <- separate(GBSXYYcoords, rn, c("pop","ind"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
 GBSXYYcoordsN$state <- GBSXYY$State
 
 GBSXYYeig <- GBSXYYpca$eig
 
 XYYPCAplot <- 
   ggplot(GBSXYYcoordsN,aes(x=Dim.2, y=Dim.1, label=pop,color=state)) + geom_text(size=10) + 
 theme_bw(base_size = 18) + guides(color = FALSE) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "PCA2 (3.17%)",y = "PCA1 (3.31%)")

 #XY
 GBSXY <- data.frame(c(GBSIDXYY[c(30:32,46:52,54:57,76:92)],GBS[c(30:32,46:52,54:57,76:92),-1]), row.names=1) #XY removed
 GBSXYpca <- PCA(GBSXY[,-(1:4)], graph = FALSE)
 fviz_pca_ind(GBSXYpca, label="", habillage=as.factor(GBSXY$State), title="", addEllipses=TRUE, ggtheme =  theme_classic())
 
 GBSXYcoord <- GBSXYpca$ind$coord
 GBSXYcoords <- setDT(data.frame(GBSXYcoord), keep.rownames = TRUE)[]
 GBSXYcoordsN <- separate(GBSXYcoords, rn, c("pop","ind"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
 GBSXYcoordsN$state <- GBSXY$State
 
 GBSXYeig <- GBSXYpca$eig

 
 XYPCAplot <- 
   ggplot(GBSXYcoordsN,aes(x=Dim.1, y=Dim.2, label=pop,color=state)) + geom_text(size=10) + 
   theme_bw(base_size = 18) + guides(color = FALSE) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "PCA1 (5.95%)",y = "PCA2 (5.03%)")
 
 multiplot(XYXYYPCAplot,XYPCAplot,XYYPCAplot,cols=3)
 
 
#,addEllipses=TRUE, ellipse.level=0.95

GBSscale <- scale(GBSs[,-(1:3)])
GBSscale[is.na(GBSscale)] <- 1 
fviz_nbclust(GBSscale, kmeans, method = "gap_stat")

GBScut <- hcut(GBSs[,-(1:3)], k = 2, stand = TRUE)
fviz_dend(GBScut, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))


#####CONSTUCT (Bradburd 2017)
library(devtools)

install_github("gbradburd/conStruct/code/conStruct",build_vignettes=TRUE)

vignette(topic="format-data",package="conStruct")
vignette(topic="run-conStruct",package="conStruct")
data(conStruct.data)

####################RNA

RNA <-  fread("allpops_aa.vclean.n012")
RNAID<-fread('name_aa.txt',header=T)
RNAs <- data.frame(c(RNAID,RNA[,-1]), row.names=1)
#RNAs <- data.frame(c(RNA[,-1]))


pwRNA<-dist.gene(RNAs,method = "pairwise")
RNAtree<-nj(as.dist(pwRNA))
RNAtree$tip.label <- RNAID$ID
plot(RNAtree,show.tip.label=TRUE,type = "unrooted", lab4ut="axial",cex=1.5,rotate.tree=-20)


RNApca <- PCåA(RNAs[,-(1:3)], graph = FALSE)
fviz_pca_ind(RNApca)
fviz_pca_ind(RNApca, label="none",habillage=as.factor(combo$Sex))
fviz_pca_ind(RNApca, label="",habillage=as.factor(RNAs$State),addEllipses=TRUE, ellipse.level=0.5)

 RNAscale <- scale(combo[,-(1:2)])
fviz_nbclust(RNAscale, kmeans, method = "gap_stat")

RNAcut <- hcut(combo[,-(1:2)], k = 3, stand = TRUE)
fviz_dend(RNAcut, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

set.seed(123)
km.res <- kmeans(RNAscale, 4, nstart = 25)

fviz_cluster(km.res, data = combo[,-(1:2)],
             ggtheme = theme_minimal(),
             main = ""
)

#sex
GBS <- fread("GBS.sex.clean.012")
GBSID<-fread('gbsnames.txt',header=T)
GBSs <- data.frame(c(GBSID,GBS[,-1]), row.names=1)

pwGBS<-dist.gene(GBSs[,-(1:3)],method = "pairwise")
tree<-nj(as.dist(pwGBS))
tree$tip.label <- GBSID$Name
plot(tree,show.tip.label=TRUE,type = "unrooted")


GBSpca <- PCA(GBSs[,-(1:3)], graph = FALSE)
fviz_pca_ind(GBSpca)
fviz_pca_ind(GBSpca, label="none")
fviz_pca_ind(GBSpca, label="GBS Population PCA",habillage=as.factor(GBSs$Pop))
fviz_pca_ind(GBSpca, label="", habillage=as.factor(GBSs$Sex),title="",addEllipses=TRUE, ellipse.level=0.95, ggtheme =  theme_classic())
