####Libraries####
setwd('~/Google Drive/Research/Data/')
options(max.print = 1000000)
library(data.table)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(sqldf)
library(GGally)
library(mclust)
library(devtools)
library(ggfortify)
library(ape)
library(FactoMineR)
library(factoextra)
library(ggridges)

####Functions####

FPKMer <- function(infile=NULL){ 
  exp <- data.frame(infile,row.names = 1)
  FPKM <-  (exp/(exp$Length/1000))/(colSums(exp, na.rm = FALSE, dims = 1)/1000000)
  FPKMt <- t(FPKM[,-c(1)])
  return(FPKMt)
}

PCAer <- function(FPKM=NULL,id=NULL){
  expPCA <- PCA(FPKM,graph=FALSE)
  print(get_eigenvalue(expPCA))
  expPCAcoord <- data.frame(expPCA$ind$coord)
  expPCAcoords <- setDT(expPCAcoord, keep.rownames = TRUE)[]
  exp_id <- fread(id)
  coords_id <- data.frame(sqldf('select exp_id.*, expPCAcoords.* 
                                from expPCAcoords left join exp_id
                                on exp_id.Individual = expPCAcoords.rn'))
  coords_id <- coords_id[,-c(6)]
  return(coords_id)
}

BIASer <- function(setOne=NULL,setTwo=NULL,locInfo=NULL,FPKM=NULL,type=NULL){
  locInfo <- fread(locInfo)
  FPKM_one <- FPKM[setOne,]
  FPKM_two <- FPKM[setTwo,]
  if (type==1) {
    FPKM_ratio <- data.frame(log((colSums(FPKM_one, na.rm = FALSE, dims = 1)/length(setOne))/
                                   (colSums(FPKM_two, na.rm = FALSE, dims = 1)/length(setTwo)),2))
    FPKM_ratio <- rename(FPKM_ratio, c("log..colSums.FPKM_one..na.rm...FALSE..dims...1..length.setOne....colSums.FPKM_two.." = "M2F"))
  }else{
    FPKM_ratio <- data.frame( ((colSums(FPKM_one, na.rm = FALSE, dims = 1)/length(setOne)) - (colSums(FPKM_two, na.rm = FALSE, dims = 1)/length(setTwo))) /
                                ((colSums(FPKM_one, na.rm = FALSE, dims = 1)/length(setOne)) + (colSums(FPKM_two, na.rm = FALSE, dims = 1)/length(setTwo)))
    )
    FPKM_ratio <- rename(FPKM_ratio, c("X..colSums.FPKM_one..na.rm...FALSE..dims...1..length.setOne....." = "M2F"))
  }
  
  FPKM_dt <- setDT(FPKM_ratio, keep.rownames = TRUE)[]
  FPKM_split <- separate(FPKM_dt, rn, c("1","locus","2"), sep = "_", remove = TRUE,
                         convert = FALSE, extra = "merge", fill = "left")
  
  FPKM_combo <- data.frame(sqldf('select FPKM_split.M2F,  locInfo.*  
                                 from locInfo left 
                                 join FPKM_split on FPKM_split.locus = locInfo.hastatulus_transcript
                                 '))
  
  FPKM_auto <- FPKM_combo[FPKM_combo$SingleAutosomal == 1,]
  FPKM_fin <- FPKM_auto[is.finite(FPKM_auto$M2F), ]
  return(FPKM_fin)
}

####Input####

allinput <- fread('exp.tsv', header=TRUE,colClasses=list(integer=2:68))

####PCA Plots####

#full set
FPKM <- FPKMer(infile=allinput)
coords_id <- PCAer(FPKM=FPKM,id='expression_id.csv')

coords_id$Population[coords_id$Population == "XYY"] <- "N"
coords_id$Population[coords_id$Population == "XY"] <- "T"

#ggpairs(coords_id[,c(4,6:10)], aes(colour = Tissue, alpha = 0.4))
#ggpairs(coords_id[,c(2,6:10)], aes(colour = Population, alpha = 0.4)) #PCA4
#ggpairs(coords_id[,c(3,6:10)], aes(colour = Sex, alpha = 0.4))

ggplot(coords_id,aes(x=Dim.1, y=Dim.2,color=Tissue,shape=Sex)) +  
  theme_bw(base_size = 18) + geom_point(size=5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5"))

ggplot(coords_id,aes(x=Dim.1, y=Dim.2,color=Tissue,shape=Sex,label=Population)) +  
  geom_text(size=5) +
  geom_point(alpha=0.5,size=10) +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5")) +
  labs(x = "PCA2 (10%)",y = "PCA3 (8.7%)") 

ggplot(coords_id,aes(x=Dim.2, y=Dim.3,color=Tissue,shape=Sex)) +  
  theme_bw(base_size = 18) + geom_point(size=5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5"))

#josh set
jFPKM <- FPKMer(infile=allinput[,c(1:2,45:68)])
coords_idj <- PCAer(FPKM=jFPKM,id='expression_id.csv')

ggplot(coords_idj,aes(x=Dim.1, y=Dim.2,color=Population,shape=Sex,label=Individual)) +  
  geom_text(size=5) + guides(shape=FALSE,color=FALSE) +
  geom_point(alpha=0.5,size=10) +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5")) + labs(title="Josh")

####BIAS ANALYSIS
#male2female Texas

jFPKM_tx <- FPKMer(infile=allinput[,c(1:2,45:56)])

jFPKM_bias <- BIASer(setOne=c(1:6),setTwo=c(7:12),locInfo='synteny_10322_hastatulus_transcripts.txt',FPKM=jFPKM_tx,type = 2)

ggplot(jFPKM_bias,aes(x=M2F)) + 
  geom_histogram(aes(y=..density..),binwidth=0.10,alpha=0.5) +
  geom_density() +
  geom_vline(xintercept = mean(mf$M2F),colour="red", linetype = "longdash") +
  scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme( panel.grid.minor = element_blank()) + 
   guides(alpha=FALSE,colour=FALSE,fill=FALSE) + labs(x="Sex-biased expression")


####Correlations####

##FMFST
fm <- fread('summarystats_FM_xy.xls')

fmpol <- melt(fm,id.vars = "locus")
fmpols <- separate(fmpol, variable, c("pop","var","cod"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
fmpols_comp <- fmpols[complete.cases(fmpols), ]
fmpols_comp$value <- as.numeric(fmpols_comp$value)

fmpols_comp$pop[fmpols_comp$pop == "pop0"] <- "All"
fmpols_comp$pop[fmpols_comp$pop == "pop1"] <- "Male"
fmpols_comp$pop[fmpols_comp$pop == "pop2"] <- "Female"

fm_fst_syn_popd <- fmpols_comp[fmpols_comp$var == 'fst' & fmpols_comp$cod == 'syn'  ,]

fm_fst_syn_split <- separate(fm_fst_syn_popd, locus, c("1","locus","2"), sep = "_", remove = TRUE,
                             convert = FALSE, extra = "merge", fill = "left")


fm_fst_syn_popd_m2fexp_xy <- data.frame(sqldf('select jFPKM_mf_bind_xy.M2F, jFPKM_mf_bind_xy.sexSpec, fm_fst_syn_split.value  
                                           from jFPKM_mf_bind_xy 
                                          left join fm_fst_syn_split on fm_fst_syn_split.locus = jFPKM_mf_bind_xy.hastatulus_transcript
                                           '))

fm_fst_bin_comp_xy <- fm_fst_syn_popd_m2fexp_xy[complete.cases(fm_fst_syn_popd_m2fexp_xy), ]


fm_fst_plot_xy <-
ggplot(fm_fst_bin_comp_xy,aes(x=M2F,y=value)) + geom_point(aes(alpha=0.5,color=sexSpec)) +
  geom_vline(xintercept = mean(fm_fst_bin_comp_xy$M2F),colour="red", linetype = "longdash") +
  geom_smooth(aes(alpha=0.7),se = FALSE,method="loess",span=0.5,size=1) + 
  theme_bw(base_size = 18)  + guides(alpha=FALSE,color=FALSE) + xlim(-3,3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="M/F Fst", x="")

ggplot(fm_fst_syn_popd_m2fexp_xy,aes(x=sexSpec,y=value)) + geom_boxplot() + labs(x="Bias Bin", y="Male to Female Fst")
fst_anova <- aov(value ~ sexSpec, data=fm_fst_syn_popd_m2fexp_xy)
summary(fst_anova) 


fm_fst_bin_xy <- summarySE(fm_fst_bin_comp_xy, measurevar="value", groupvars=c("sexSpec"))


fm_fst_bin_plot_xy <-
ggplot(fm_fst_bin_xy, aes(x=sexSpec, y=value, fill=sexSpec)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + 
  labs(x="Expression Bias", y="M/F Fst") +
  scale_x_discrete(limits=c("extFem","intFem","None","intMale","extMale"))


multiplot(fm_fst_plot_xy,fm_fst_bin_plot_xy)

#FMdxy
fm_dxy_syn_popd <- fmpols_comp[fmpols_comp$var == 'Dxy' & fmpols_comp$cod == 'syn'  ,]

fm_dxy_syn_split <- separate(fm_dxy_syn_popd, locus, c("1","locus","2"), sep = "_", remove = TRUE,
                             convert = FALSE, extra = "merge", fill = "left")

dxy_syn_popd_bias <- data.frame(sqldf('select jFPKM_mf_bind_xy.M2F, jFPKM_mf_bind_xy.sexSpec, fm_dxy_syn_split.value  
                                           from jFPKM_mf_bind_xy 
                                            left join fm_dxy_syn_split on fm_dxy_syn_split.locus = jFPKM_mf_bind_xy.hastatulus_transcript
                                            '))

dxy_syn_popd_bias_comp <- dxy_syn_popd_bias[complete.cases(dxy_syn_popd_bias), ]


fm_dxy_plot <- 
ggplot(dxy_syn_popd_bias_comp,aes(x=M2F,y=value)) + geom_point(aes(alpha=0.5,color=sexSpec)) +
  geom_vline(xintercept = mean(dxy_syn_popd_bias_comp$M2F),colour="red", linetype = "longdash") +
  geom_smooth(aes(alpha=0.7),se = FALSE,method="loess",span=0.5,size=2) + 
  theme_bw(base_size = 18)  + guides(alpha=FALSE,color=FALSE) + xlim(-3,3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="M/F Dxy", x="")


dxy_anova <- aov(value ~ sexSpec, data=dxy_syn_popd_bias_comp)
summary(dxy_anova) 


dxy_bin <- summarySE(dxy_syn_popd_bias_comp, measurevar="value", groupvars=c("sexSpec"))

fm_dxy_bin_plot <- 
ggplot(dxy_bin, aes(x=sexSpec, y=value, fill=sexSpec)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Expression Bias", y="M/F Dxy") +
  scale_x_discrete(limits=c("extFem","intFem","None","intMale","extMale"))


multiplot(fm_dxy_plot,fm_dxy_bin_plot)
multiplot(fm_fst_plot_xy,fm_fst_bin_plot_xy,fm_dxy_plot,fm_dxy_bin_plot,cols=2)


#TajD
tajd_syn_popd <- pols_comp[pols_comp$var == 'TajD' & pols_comp$cod == 'syn' ,]

tajd_syn_popd_bias <- data.frame(sqldf('select jFPKM_mf_bind_xy.M2F, jFPKM_mf_bind_xy.sexSpec, tajd_syn_popd.value, tajd_syn_popd.pop  
                                      from jFPKM_mf_bind_xy 
                                      left join tajd_syn_popd on  tajd_syn_popd.locus = jFPKM_mf_bind_xy.hastatulus_transcript
                                       '))


tajd_bin <- summarySE( tajd_syn_popd_bias_auto, measurevar="value", groupvars=c("sexSpec","pop"))

tajd_bin_sh <- tajd_bin[tajd_bin$pop != "All",]

sexTajDplot <- 
ggplot(tajd_bin_sh, aes(x=sexSpec, y=value, fill=pop)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Expression Bias", y="TajD") + 
  scale_fill_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5")) +
  scale_x_discrete(limits=c("extFem","intFem","none","intMale","extMale"))

tajd_syn_popd_bias_auto <- tajd_syn_popd_bias_auto[tajd_syn_popd_bias_auto$pop != "All",]

tajd_syn_popd_bias_anova <- aov(value ~ sexSpec*pop, data=tajd_syn_popd_bias_auto)

summary(tajd_syn_popd_bias_anova) # display Type I ANOVA table



####Pollen2Flower####
gFPKM <- t(FPKM[,c(2:5,10:17)])
gFPKM_pol <- gFPKM[c(1:4),]
gFPKM_flow <- gFPKM[-c(1:4),]

#study pca
expPCAg <- PCA(gFPKM,graph=FALSE)
expPCAcoordg <- data.frame(expPCAg$ind$coord)
expPCAcoordsg <- setDT(expPCAcoordg, keep.rownames = TRUE)[]
coords_idg <- data.frame(sqldf('select exp_id.*, expPCAcoordsg.* 
                               from expPCAcoordsg left join exp_id
                               on exp_id.Individual = expPCAcoordsg.rn'))

coords_idg <- coords_idg[,-c(6)]

#ggpairs(coords_idg[,c(4,6:10)], aes(colour = Tissue, alpha = 0.4))

#josh <- 
ggplot(coords_idg,aes(x=Dim.1, y=Dim.2,color=Population,shape=Sex,label=Individual)) +  
  geom_text(size=5) + guides(shape=FALSE,color=FALSE) +
  geom_point(alpha=0.5,size=10) +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5")) + labs(title="Josh")


#bias
gFPKM_pf <- data.frame(log((colSums(gFPKM_pol, na.rm = FALSE, dims = 1)/4)/
                             (colSums(gFPKM_flow, na.rm = FALSE, dims = 1)/4),2))
gFPKM_pf <- rename(gFPKM_pf, c("log..colSums.gFPKM_pol..na.rm...FALSE..dims...1..4...colSums.gFPKM_flow.." = "P2F"))
gFPKM_pf_dt <- setDT(gFPKM_pf, keep.rownames = TRUE)[]


#m2f_plot <-
ggplot(gFPKM_pf,aes(x=P2F)) + 
  geom_histogram(aes(y=..density..),binwidth=0.50,alpha=0.5) +
  geom_density() +
  scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme( panel.grid.minor = element_blank()) + 
  guides(alpha=FALSE,colour=FALSE,fill=FALSE) + labs(x="log Male to Female expression")

#mclusters
gFPKM_pf_comp <- gFPKM_pf_dt[complete.cases(gFPKM_pf_dt) ]
gFPKM_pf_compl <- gFPKM_pf_comp[!is.infinite(gFPKM_pf_comp$P2F) ]

gFPKM_pf_clust <- densityMclust(gFPKM_pf_compl$P2F,G=3)
gFPKM_pf_clust_m <- Mclust(gFPKM_pf_compl$P2F,G=3)

#par(mfrow = c(2,1))
plot(gFPKM_pf_clust,what="density",data=gFPKM_pf_compl$P2F,breaks=50)
#plot(jFPKM_mf_clust_m,what="classification")
par(mfrow = c(1,1))

class1 <- gsub("1","Flower",gFPKM_pf_clust$classification)
class2 <- gsub("2","None",class1)
class3 <- gsub("3","Pollen",class2)


gFPKM_pf_bind <- cbind(gFPKM_pf_compl,as.factor(class3))

gFPKM_pf_split <- separate(gFPKM_pf_bind, rn, c("1","locus","2"), sep = "_", remove = TRUE,
                           convert = FALSE, extra = "merge", fill = "left")


ggplot(gFPKM_pf_bind,aes(x=V2,y=P2F)) + geom_boxplot() 


##correlations
##tajD
tajd_syn_popd_bias_pf <- data.frame(sqldf('select gFPKM_pf_split.locus, gFPKM_pf_split.P2F, gFPKM_pf_split.V2, tajd_syn_popd.value, tajd_syn_popd.pop  
                                       from gFPKM_pf_split left join tajd_syn_popd
                                       on  tajd_syn_popd.locus = gFPKM_pf_split.locus'))

tajd_syn_popd_bias_pf_comp <- tajd_syn_popd_bias_pf[complete.cases(tajd_syn_popd_bias_pf), ]

tajd_bin_pf <- summarySE( tajd_syn_popd_bias_pf_comp, measurevar="value", groupvars=c("V2","pop"))

tajd_bin_pf_sh <- tajd_bin_pf[tajd_bin_pf$pop != "All",]

pollenTajDplot <-
ggplot(tajd_bin_pf_sh, aes(x=V2, y=value, fill=pop)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Expression Bias", y="") + 
  scale_fill_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5"))

tajd_syn_popd_bias_pf_comp <- tajd_syn_popd_bias_pf_comp[tajd_syn_popd_bias_pf_comp$pop != "All",]

tajd_syn_popd_bias_pf_anova <- aov(value ~ V2*pop, data=tajd_syn_popd_bias_pf_comp)

summary(tajd_syn_popd_bias_pf_anova)


multiplot(sexTajDplot,pollenTajDplot, cols=2)

##dnds

D_syn_popd <- pols_comp[pols_comp$var == 'Dxy' & pols_comp$cod == 'syn' ,]
D_rep_popd <- pols_comp[pols_comp$var == 'Dxy' & pols_comp$cod == 'rep' ,]
pi_syn_popd <- pols_comp[pols_comp$var == 'pi' & pols_comp$cod == 'syn' ,]
pi_rep_popd <- pols_comp[pols_comp$var == 'pi' & pols_comp$cod == 'rep' ,]




mkd_popd_bias <- data.frame(sqldf('select jFPKM_mf_bind_xy.hastatulus_transcript, jFPKM_mf_bind_xy.M2F, jFPKM_mf_bind_xy.sexSpec,
                           D_syn_popd.value, 
                          D_rep_popd.value, 
                          pi_syn_popd.value, pi_syn_popd.pop,
                          pi_rep_popd.value, pi_rep_popd.pop
                                       from jFPKM_mf_bind_xy 
                                      left join D_syn_popd
                                       on  jFPKM_mf_bind_xy.hastatulus_transcript = D_syn_popd.locus
                                         left join  D_rep_popd
                                       on  jFPKM_mf_bind_xy.hastatulus_transcript =  D_rep_popd.locus
                                        left join  pi_syn_popd
                                       on  jFPKM_mf_bind_xy.hastatulus_transcript =  pi_syn_popd.locus
                                         left join  pi_rep_popd
                                       on  jFPKM_mf_bind_xy.hastatulus_transcript =  pi_rep_popd.locus
                                        
                                  '))

mkd_popd_bias$dnds <- mkd_popd_bias$value.1/mkd_popd_bias$value
mkd_popd_bias$pnps <- mkd_popd_bias$value.3/mkd_popd_bias$value.2


mkd_popd_bias <- mkd_popd_bias[c(3,4:11)]

dnds_bin_comp <- mkd_popd_bias[complete.cases(mkd_popd_bias), ]
pnps_bin_comp <- dnds_bin_comp[!is.infinite(dnds_bin_comp$pnps), ]
pnps_bin_popd <- pnps_bin_comp[pnps_bin_comp$pop.1 == pnps_bin_comp$pop,]


dnds_bin_anova <- aov(value ~ sexSpec, data=dnds_bin_comp)
summary(dnds_bin_anova)

pnps_bin_anova <- aov(value ~ sexSpec*pop, data=pnps_bin_popd)
summary(pnps_bin_anova)


dnds_bin <- summarySE( dnds_bin_comp, measurevar="dnds", groupvars=c("sexSpec"))
pnps_bin <- summarySE( pnps_bin_popd, measurevar="pnps", groupvars=c("sexSpec","pop"))

dnds_plot <-
ggplot(dnds_bin, aes(x=sexSpec, y=dnds,fill=sexSpec)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=dnds-se, ymax=dnds+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Expression Bias", y="Dn/Ds") +
  scale_x_discrete(limits=c("extFem","intFem","None","intMale","extMale"))

title <- expression(paste(pi, ""[n],"/", pi, ""[s]))

pnps_bin_sh <- pnps_bin[pnps_bin$pop != "All",]

pnps_plot <-
ggplot(pnps_bin_sh, aes(x=sexSpec, y=pnps,fill=pop)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=pnps-se, ymax=pnps+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Expression Bias", y=title) + 
  scale_fill_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5"))+
  scale_x_discrete(limits=c("extFem","intFem","None","intMale","extMale"))


multiplot(dnds_plot,pnps_plot)

#PollenFlower

mkd_popd_PFbias <- data.frame(sqldf('select gFPKM_pf_split.locus, gFPKM_pf_split.P2F, gFPKM_pf_split.V2, locInfo.*,
                           D_syn_popd.value, 
                                  D_rep_popd.value, 
                                  pi_syn_popd.value, pi_syn_popd.pop,
                                  pi_rep_popd.value, pi_rep_popd.pop
                                  from gFPKM_pf_split 
                                  left join D_syn_popd
                                  on  gFPKM_pf_split.locus = D_syn_popd.locus
                                  left join  D_rep_popd
                                  on  gFPKM_pf_split.locus =  D_rep_popd.locus
                                  left join  pi_syn_popd
                                  on  gFPKM_pf_split.locus =  pi_syn_popd.locus
                                  left join  pi_rep_popd
                                  on  gFPKM_pf_split.locus =  pi_rep_popd.locus
                                  left join locInfo
                                  on gFPKM_pf_split.locus = locInfo.hastatulus_transcript
                                    
                                  '))

mkd_popd_PFbias$dnds <- mkd_popd_PFbias$value.1/mkd_popd_PFbias$value
mkd_popd_PFbias$pnps <- mkd_popd_PFbias$value.3/mkd_popd_PFbias$value.2

mkd_popd_PFbias_auto <- mkd_popd_PFbias[mkd_popd_PFbias$TXjAuto == 1,]
mkd_popd_PFbias_auto <- mkd_popd_PFbias_auto[c(3,13:20)]

dnds_pfbin_comp <- mkd_popd_PFbias_auto[complete.cases(mkd_popd_PFbias_auto), ]
pnps_pfbin_comp <- dnds_pfbin_comp[!is.infinite(dnds_pfbin_comp$pnps), ]
pnps_pfbin_popd <- pnps_pfbin_comp[pnps_pfbin_comp$pop.1 == pnps_pfbin_comp$pop,]

dnds_pfbin_anova <- aov(value ~ V2, data=dnds_pfbin_comp)
summary(dnds_pfbin_anova)

pnps_pfbin_anova <- aov(value ~ V2*pop, data=pnps_pfbin_popd)
summary(pnps_pfbin_anova)

dnds_pfbin <- summarySE( dnds_pfbin_comp, measurevar="dnds", groupvars=c("V2"))
pnps_pfbin <- summarySE( pnps_pfbin_popd, measurevar="pnps", groupvars=c("V2","pop"))


dnds_pfplot <-
  ggplot(dnds_pfbin, aes(x=V2, y=dnds,fill=V2)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=dnds-se, ymax=dnds+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="", y="") 

pnps_pfbin <- pnps_pfbin[pnps_pfbin$pop != "All",  ]

pnps_pfplot <-
  ggplot(pnps_pfbin, aes(x=V2, y=pnps,fill=pop)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=pnps-se, ymax=pnps+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Bias", y="") + 
  scale_fill_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5"))


multiplot(dnds_pfplot,pnps_pfplot)
multiplot(dnds_plot,pnps_plot,dnds_pfplot,pnps_pfplot,cols=2)

####LD#####

ld_table <- 
  data.frame(cbind(bias = c("None","Low","High"),mean = c(0.09090,0.07639,0.07088),se = c(3.106e-06,6.743e-06,0.0002116)))

ld_table$mean <- as.numeric(as.character(ld_table$mean))
ld_table$se <- as.numeric(as.character(ld_table$se))


ggplot(ld_table, aes(x=bias, y=mean,fill=bias)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x="Sex Expresion Bias", y="LD") #+ 
  scale_fill_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5"))

  
  
####expMK = deltaX####

#r = divergence/variance

jFPKM <- t(FPKM[,c(44:67)])
jFPKM_xy <- jFPKM[c(1:12),]
jFPKM_xyy <- jFPKM[-c(1:12),]

#xy
jFPKM_xy_t <- data.frame(t(jFPKM_xy))
jFPKM_xy_dt <- setDT(jFPKM_xy_t, keep.rownames = TRUE)[]

jFPKM_xy_melt <- melt(jFPKM_xy_dt,id.vars = "rn")
jFPKM_xy_summary <- summarySE( jFPKM_xy_melt, measurevar="value", groupvars=c("rn"))
jFPKM_xy_summary$vari <- jFPKM_xy_summary$sd / jFPKM_xy_summary$value

ggplot(jFPKM_xy_summary,aes(x=vari)) + geom_density() 

jFPKM_xy_comp <- jFPKM_xy_summary[complete.cases(jFPKM_xy_summary), ]
sum(jFPKM_xy_comp$vari) / length(jFPKM_xy_comp$vari)
#0.973

#xyy
jFPKM_xyy_t <- data.frame(t(jFPKM_xyy))
jFPKM_xyy_dt <- setDT(jFPKM_xyy_t, keep.rownames = TRUE)[]

jFPKM_xyy_melt <- melt(jFPKM_xyy_dt,id.vars = "rn")
jFPKM_xyy_summary <- summarySE( jFPKM_xyy_melt, measurevar="value", groupvars=c("rn"))
jFPKM_xyy_summary$vari <- jFPKM_xyy_summary$sd / jFPKM_xyy_summary$value

ggplot(jFPKM_xyy_summary,aes(x=vari)) + geom_density() 

jFPKM_xyy_comp <- jFPKM_xyy_summary[complete.cases(jFPKM_xyy_summary), ]
sum(jFPKM_xyy_comp$vari) / length(jFPKM_xyy_comp$vari)
#0.865

jFPKM_variances <- cbind(XYY=jFPKM_xyy_summary$vari,XY=jFPKM_xy_summary$vari)

jFPKM_variances_melt <- melt(jFPKM_variances)

jFPKM_variances_complete <- jFPKM_variances_melt[complete.cases(jFPKM_variances_melt),]

ggplot(jFPKM_variances_complete, aes(x = value, y = Var2)) + geom_density_ridges()


for (locus in c(1:10000)){

  
  nXY = length(XYsamples)
  nXYY = length(XYYsamples) 
  
  for (indXY in c(1:nXY)){
  
    
    var()
  
    
  }
  
  for (indXY in c(1:nXY)){
  
    divXtemp = log( FPKM$indXY / FPKM$indXYY   
      ,2)
  
    divXtot = divXtot + divXtemp
  }
  divX = divXtot / n 
  
  
  rXxyy = divX / varXxyy
  
}




####StudyBias#### 
fFPKM <- t(FPKM[,c(22:43,68)])

expPCAf <- PCA(fFPKM,graph=FALSE)
expPCAcoordf <- data.frame(expPCAf$ind$coord)
expPCAcoordsf <- setDT(expPCAcoordf, keep.rownames = TRUE)[]
coords_idf <- data.frame(sqldf('select exp_id.*, expPCAcoordsf.* 
                               from exp_id left join expPCAcoordsf 
                               on exp_id.Individual = expPCAcoordsf.rn'))

coords_idf <- coords_idf[,-c(6)]

felix <-
  ggplot(coords_idf,aes(x=-(Dim.2), y=Dim.1,color=Population,shape=Sex,label=Individual)) + 
  geom_text(size=5) + guides(shape=FALSE,color=FALSE) +
  geom_point(alpha=0.5,size=10) +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values=c( "#ffd700", "#fa8775",   "#cd34b5")) + labs(title="Felix")

multiplot(felix,georg,josh,cols=3)



