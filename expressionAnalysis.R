#### TOP ####
setwd('~/Google Drive/Research/Data/RNAseq')
library(data.table)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(sqldf)
library(rapport)
library(ggridges)
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("DESeq2"))
#biocLite(c("apeglm"))
library(DESeq2)
library(apeglm)
library(FactoMineR)
library(factoextra)
options(max.print = 1000000)


####HomeMadeFunctions####
statsExp <- function(set=NULL, outgroup="rothschildianus", popStr="pop", chrom=NULL, highbias=NULL, lowbias=NULL) {
  
  filename_win <- paste(set,"_",outgroup,"_summarystats_",popStr,"_",chrom,".txt",sep="")
  filename_bw <- paste(set,"_",outgroup,"_interpop_",popStr,"_",chrom,".txt",sep="")
  inter <- fread(filename_bw)
  within <- fread(filename_win)
  inter <- separate(inter, locus, c("locus","file"), 
                    sep = ".fasta", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
  within <- separate(within, locus, c("locus","file"), 
                     sep = ".fasta", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
  
  interExp <- data.frame(sqldf('select resdf.*, inter.*, within.* from resdf 
                             left join inter on resdf.loci = inter.locus 
                               left join within on resdf.loci = within.locus'))
  
  interExp$sexBias[interExp$log2FoldChange > 0] <- highbias
  interExp$sexBias[interExp$log2FoldChange < 0] <- lowbias
  interExp$sexBias[interExp$pvalue >  0.05] <- "Unbiased"
  
  interExp <- interExp[!is.na(interExp$locus),]
}


####import####

base <- fread('allLoci.list',header = FALSE)
inds <- fread('rna.sql.list',header = FALSE)
inds <- inds$V1
inds <- c(inds,"SC1pollen","SC2pollen","TX1Bpollen","TX2Bpollen")
inds <- c(inds,"9_SCFemFlower","10_SCFemFlower","11_SCFemFlower","12_SCFemFlower",
          "13_SCMaleFlower","14_SCMaleFlower","15_SCMaleFlower","16_SCMaleFlower",
          "17_TXMaleFlower","18_TXMaleFlower","19_TXMaleFlower","20_TXMaleFlower",
          "21_TXFemFlower","22_TXFemFlower","23_TXFemFlower","24_TXFemFlower")

#MERGE#
loops <- 1
reads <- base

for (ind in inds){
  cat(ind,"\n")
  ind <- toString(ind)
  fileName <- paste(ind,".reads.txt",sep="")
  file <- fread(fileName,header = FALSE)
  names(file) <- c("loci","length",ind,"unmapped")
  reads <- data.frame(sqldf('select reads.*, file.* from reads left join file on reads.V1 = file.loci'))
  columns <- c(1:(loops+1),(loops+4))
  if (loops == 1){
    columns <- c(1,3,4)
  }
  loops <- loops + 1
  reads <- reads[,columns]
}

####Orthologs####

roth <- fread('rothschildianus.orthologs',header=FALSE)
buc <- fread('bucephalophorus.orthologs',header=FALSE)

species <- data.frame(sqldf('select base.*, roth.*, buc.* from base left join roth on base.V1 = roth.V1 left join buc on base.V1 = buc.V1'))
species <- species[,c(1,3,5)]
names(species) <- c("hast","roth","buc")

rothInds <- c("roth_cross_1_Female_mother","roth_cross_1_Male_father","roth_cross_2_Female_mother","roth_cross_2_Male_father",
"Rothpollen","1_rothFemFlower","2_rothFemFlower","3_rothFemFlower","4_rothFemFlower",
"5_rothFemFlower","6_rothMaleFlower","7_rothMaleFlower","8_rothMaleFlower")

loops <- 1
convert <- species
for (ind in rothInds){
  #remind me where I am
  cat(ind,"\n")
  ind <- toString(ind)
  #import file
  fileName <- paste(ind,".reads.txt",sep="")
  file <- fread(fileName,header = FALSE)
  names(file) <- c("loci","length",ind,"unmapped")
  #add on reads
  convert <- data.frame(sqldf('select convert.*, file.* from convert left join file on convert.roth = file.loci'))
  columns <- c(1:(loops+2),(loops+5))
  if (loops == 1){
    columns <- c(1,2,3,6)
  }
  loops <- loops + 1
  convert <- convert[,columns]
}

bucInds <- c("6.RB1leaf")
for (ind in bucInds){
  #remind me where I am
  cat(ind,"\n")
  ind <- toString(ind)
  #import file
  fileName <- paste(ind,".reads.txt",sep="")
  file <- fread(fileName,header = FALSE)
  names(file) <- c("loci","length",ind,"unmapped")
  #add on reads
  convert <- data.frame(sqldf('select convert.*, file.* from convert left join file on convert.buc = file.loci'))
  columns <- c(1:(loops+2),(loops+5))
  if (loops == 1){
    columns <- c(1,2,3,6)
  }
  loops <- loops + 1
  convert <- convert[,columns]
}

allExp <- data.frame(sqldf('select  reads.*, convert.* from reads left join convert on reads.V1 = convert.hast '))
allExp <- allExp[,-c(70:72)]
allExp <- allExp[complete.cases(allExp),] 

pop <- c("T","T","T","T","T","T","T","T","T","T","T","T","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","T","N","N","N","N","N","T","T","N","N","N","T","T","T","T","T","N","N","T","T","N","N","N","N","N","N","N","N","T","T","T","T","T","T","T","T","R","R","R","R","R","R","R","R","R","R","R","R","R","B")
sex <- c("M","M","M","F","F","F","F","F","F","M","M","M","F","F","F","F","F","F","M","M","M","M","M","M","M","M","I","M","M","M","M","M","M","I","I","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","F","F","F","F","M","M","M","M","M","M","M","M","F","F","F","F","F","M","F","M","M","F","F","F","F","F","M","M","M","H") ##Add weird males
origin <- c("J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","D","D","D","D","G","G","G","G","G","G","G","G","G","J")
tissue <- c("L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","P","P","P","P","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","L","L","L","L","P","F","F","F","F","F","F","F","F","L")

####PCA####

readsCounts <- allExp[,-c(1)]
FPKM <-  data.frame(t( t(readsCounts/  (readsCounts$length/1000)  )/(colSums(readsCounts, na.rm = FALSE, dims = 1)/1000000) ) )
FPKM <- t(FPKM[,-c(1)])
ExpPCA <- PCA(FPKM, graph = FALSE)

allinds <- c(inds,rothInds,bucInds)
countInds<-c(1:length(allinds))
infor <- cbind(allinds,countInds,sex,pop,origin,tissue)

ExpPCAcoord <- ExpPCA$ind$coord
bound <- data.frame(cbind(infor,ExpPCAcoord))

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

bound$Dim1 <- as.numeric.factor(bound$Dim.1)
bound$Dim2 <- as.numeric.factor(bound$Dim.2)

ExpPCAeig <- ExpPCA$eig

tPCA <- ggplot(bound,aes(x=Dim1, y=Dim2, label=countInds,color=tissue)) + geom_text(size=10) + 
  theme_bw(base_size = 18) + #guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title="tissue")

sPCA <- ggplot(bound,aes(x=Dim1, y=Dim2, label=countInds,color=sex)) + geom_text(size=10) + 
  theme_bw(base_size = 18) + #guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title="sex")

pPCA <- ggplot(bound,aes(x=Dim1, y=Dim2, label=countInds,color=pop)) + geom_text(size=10) + 
  theme_bw(base_size = 18) + #guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title="pop")

oPCA <- ggplot(bound,aes(x=Dim1, y=Dim2, label=countInds,color=origin)) + geom_text(size=10) + 
  theme_bw(base_size = 18) + #guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title="origin")

multiplot(tPCA,sPCA,pPCA,oPCA,cols=2)

####DESEQ####
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#Read Intake#
treads <- allExp[,-c(1,2)]
rownames(treads) <- allExp$V1
countdata <- data.matrix(treads,rownames.force=TRUE) 
colData <- cbind(sex,pop,origin,tissue)
rownames(colData) <- colnames(allExp[,-c(1,2)])

##Define bias analysis##
#dds <- DESeqDataSetFromMatrix(countData = treads, colData = colData, design = ~ sex)
#dds <- DESeqDataSetFromMatrix(countData = treads, colData = colData, design = ~  pop)
dds <- DESeqDataSetFromMatrix(countData = treads, colData = colData, design = ~ sex + tissue )

dds <- dds[ rowSums(counts(dds)) > 1, ] #filter for rows with info
dss <- DESeq(dds)

#res <- results(dss, contrast=c("pop","T","N"),  alpha = 0.05  )
#res <- results(dss, contrast=c("sex","M","F"),  alpha = 0.05  )
res <- results(dss, contrast=c("tissue","L","P"),  alpha = 0.05  )
resdf <- data.frame(res)
resdf$loci <- row.names(resdf)
interExp <- resdf

#plotMA(res, ylim=c(-2,2))

####listExport####

##export lists of biased genes##
#pollenList <- interExp$loci[interExp$sexBias == "Pollen"]
#write(pollenList, file = "pollen.list",
#      append = FALSE, sep = "\n")

#leafList <- interExp$loci[interExp$sexBias == "Leaf"]
#write(leafList, file = "leaf.list",
#      append = FALSE, sep = "\n")

#LPUnbiasedList <- interExp$loci[interExp$sexBias == "LPUnbiased"]
#write(LPUnbiasedList, file = "LPUnbiased.list",
#      append = FALSE, sep = "\n")

####regions####

alist <- fread('A.list',header=FALSE)
xylist <- fread('XY.list',header=FALSE)
hlist <- fread('H.list',header=FALSE)
nlist <- fread('N.list',header=FALSE)

regions <- data.frame(sqldf('select base.*, alist.*, xylist.*, hlist.* , nlist.* 
                            from base 
                            left join alist on base.V1 = alist.V1 
                            left join xylist on base.V1 = xylist.V1
                            left join hlist on base.V1 = hlist.V1 
                            left join nlist on base.V1 = nlist.V1
                            '))
names(regions) <- c("Loci","A","XY","H","N")

regions$regions[ !is.na(regions$A) == TRUE ] <- "Autosomal"
regions$regions[ !is.na(regions$XY) == TRUE ] <- "XY"
regions$regions[ !is.na(regions$H) == TRUE ] <- "Hemizygous"
regions$regions[ !is.na(regions$N) == TRUE ] <- "NeoXY"
regions <- regions[,c(1,6)]

regionsExp <- data.frame(sqldf('select resdf.*, regions.* from resdf left join regions on resdf.loci = regions.Loci'))
regionsExp <- regionsExp[,-8]
#regionsExp <- regionsExp[regionsExp$pvalue <  0.05,]
#regionsExp <- regionsExp[!is.na(regionsExp$pvalue), ]


ggplot(regionsExp, aes(x = log2FoldChange, y = regions, fill=regions, alpha=0.5)) + 
  geom_density_ridges() + guides(alpha=FALSE,fill=FALSE) +
  theme_bw()  + theme_bw(base_size = 30) + 
 # labs(x = "log2 Male Expression Bias", y="") +
  labs(x = "log2 Leaf Expression Bias", y="") +
scale_fill_manual(values=c( 
  '#00ADEF', #Blue
  '#FFF100',  #yellow
  '#00A550', #green
 # '#1B75BB', #purple-y
  '#8B4BD8'
)) + 
  xlim(-3,3)

####regionStacks####

stacks <- data.frame(
cbind(
  rbind(
    length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]) /  nrow(xylist), #/ (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ])),
    length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]) /  nrow(xylist)#/ (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]))
  ),
  rbind(
    length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ] ) /  nrow(alist), # / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ])),
    length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ] ) /  nrow(alist)#/ (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ]))
  ),
  rbind(
    length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ] )/  nrow(hlist), #/ (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ])),
    length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ] )/  nrow(hlist) #/ (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ]))
  ),
  rbind(
    length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ] )/  nrow(nlist), # / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ])),
    length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ] )/  nrow(nlist) #/ (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ]))
  )
)
#,row.names = c("Female","Male"))
,row.names = c("Pollen","Leaf"))
#,row.names = c("Texas","NC"))


names(stacks) <- c("XY","Autosomal","Hemizygous","NeoXY")
tstacks <- data.frame(t(stacks))
tstacks$region <- rownames(tstacks)
stack_melt <- melt(tstacks,id.vars = "region")

ggplot(stack_melt, aes(x = region, y = value,fill=variable)) + 
  geom_bar(stat="identity") +
  # facet_grid(. ~ regions) +
  labs(x = "", y="Proportion with bias expression",fill="") +
  theme_bw()  + theme_bw(base_size = 30) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    #'#FFF100',  #yellow
    #'#00A550', #green
    # '#1B75BB', #purple-y
    '#8B4BD8'
  )) #+
  ylim(-.5,0.5)

####Stat Correlates####

#Tajima's D#

ggplot(interExp, aes(x = log2FoldChange, y = pop0_tajD_syn, color=pvalue)) + 
  geom_point() + 
  stat_smooth(method = "loess") +
 # facet_grid(. ~ regions) +
  labs(x = "log2 Expression Bias", y="TajD_syn") +
  theme_bw()  + theme_bw(base_size = 30) + 
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    # '#1B75BB', #purple-y
    '#8B4BD8'
  ))

title_tds <- expression(paste(TajD, ""[syn]))

ggplot(interExp, aes(y=pop0_tajD_syn, x=sexBias)) + 
  geom_violin() + geom_boxplot(width=0.1) +
  theme_bw()  + theme_bw(base_size = 30) +
  #facet_grid(. ~ regions) +
  labs(x = "Expression Bias", y=title_tds) #+
  scale_x_discrete(limits=c("Male","Unbiased","Female")) 

t.test(interExp$pop0_tajD_syn[interExp$sexBias == "Male"],interExp$pop0_tajD_syn[interExp$sexBias == "Female"])

TajD_LP_anova <- aov(pop0_tajD_syn ~ sexBias, data=interExp)
summary(TajD_LP_anova) 

#pi#

title_pisyn <- expression(paste(pi, ""[syn]))

ggplot(interExp, aes(y=pop0_pi_syn, x=sexBias)) + 
  geom_violin() + geom_boxplot(width=0.1) +
  theme_bw()  + theme_bw(base_size = 30) +
  #facet_grid(. ~ regions) +
  labs(x = "Expression Bias", y=title_pisyn) #+
  scale_x_discrete(limits=c("Male","Unbiased","Female")) +
  ylim(0,0.01)

t.test(interExp$pop0_pi_syn[interExp$sexBias == "Male"],interExp$pop0_pi_syn[interExp$sexBias == "Female"])

pi_MF_anova <- aov(pop0_pi_syn ~ sexBias, data=interExp)
summary(pi_MF_anova) 

interExpBiased <- interExp[interExp$pvalue < 0.05,]
interExpBiased$log2FoldChangeAbs <- abs(interExpBiased$log2FoldChange)

interExpBiased_comp <- interExpBiased[!is.na(interExpBiased$sexBias), ]

ggplot(interExpBiased_comp, aes(y=pop0_pi_syn, x=log2FoldChangeAbs, color=sexBias)) + 
  geom_point(alpha=0.1) + guides(alpha=FALSE) +
  stat_smooth(method = "lm",se=FALSE,size=2) +
  theme_bw()  + theme_bw(base_size = 30) +
  labs(x = "Expression Bias", y=title_pisyn, color="Sex") +
  ylim(0,0.01)

pi_male_fit <- lm(pop0_pi_syn ~ log2FoldChangeAbs , data=interExpBiased_comp[interExpBiased_comp$sexBias == "Male",])
summary(pi_male_fit) 

pi_female_fit <- lm(pop0_pi_syn ~ log2FoldChangeAbs , data=interExpBiased_comp[interExpBiased_comp$sexBias == "Female",])
summary(pi_female_fit) 



##dnds##

ggplot(interExp, aes(y=pop0_dnds_NA, x=log2FoldChange, color=sexBias)) + 
  geom_point() + guides(alpha=FALSE) +
  #stat_smooth(method = "lm",se=FALSE,size=2) +
  theme_bw()  + theme_bw(base_size = 30) +
  labs(x = "Expression Bias", y="dnds", color="") 

#A_dnds <-
ggplot(interExp, aes(y=pop0_dnds_NA, x=sexBias)) + 
  #geom_violin() + 
  geom_boxplot(width=0.1) +
  #geom_point(aes())+
  stat_smooth(method = "lm",se=FALSE) + guides(color=FALSE) +
  theme_bw()  + theme_bw(base_size = 30) + #xlim(0,10) + #ylim(0,5) + 
  #facet_grid(. ~ regions) + 
  labs(x = "log Expression difference", y="dn/ds", color="Bias") +
  scale_x_discrete(limits=c("Leaf","Unbiased","Pollen")) #+
  #scale_x_discrete(limits=c("Male","Unbiased","Female"))
  annotate(geom="text", x = "Pollen", y = 2, label = "b", parse = TRUE, size=10) +
  annotate(geom="text", x = "Unbiased", y = 2, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "Leaf", y = 2, label = "ab", parse = TRUE, size=10) 
  

t.test(interExp$pop0_dnds_NA[interExp$sexBias == "Leaf"],interExp$pop0_dnds_NA[interExp$sexBias == "Pollen"])
t.test(interExp$pop0_dnds_NA[interExp$sexBias == "Leaf"],interExp$pop0_dnds_NA[interExp$sexBias == "Unbiased"])
t.test(interExp$pop0_dnds_NA[interExp$sexBias == "Unbiased"],interExp$pop0_dnds_NA[interExp$sexBias == "Pollen"])


mean(!is.na(interExp$pop0_dnds_NA[interExp$sexBias == "Leaf"]))
mean(!is.na(interExp$pop0_dnds_NA[interExp$sexBias == "Pollen"]))

####mk####

title_mk <- expression(paste("dn/ds/",pi, ""[n],"/",pi, ""[s]))

interExp$log2FoldChangeAbs <- abs(interExp$log2FoldChange)
interExp$log2FoldChangeAbsRound <- as.factor(round(interExp$log2FoldChangeAbs,0))

statsExp(set="m",chrom="Y",highbias="Leaf",lowbias="Pollen")

#X_mk <-
#Y_mk <- 
ggplot(interExp, aes(y=pop0_mk_NA,  x=sexBias,color=sexBias)) + 
  #geom_violin() + 
  geom_boxplot(width=0.1) + 
  #geom_point(aes())+
  stat_smooth(method = "lm",se=FALSE) + guides(color=FALSE) +
  theme_bw()  + theme_bw(base_size = 30) +  ylim(0,4) + #xlim(0,10)  + 
  #facet_grid(. ~ regions) + 
  labs(x = "", y=title_mk, color="Bias",title="X") +
  scale_x_discrete(limits=c("Leaf","Unbiased","Pollen")) +
  #scale_x_discrete(limits=c("Male","Unbiased","Female"))
  annotate(geom="text", x = "Pollen", y = 2, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "Unbiased", y = 2, label = "ab", parse = TRUE, size=10) +
  annotate(geom="text", x = "Leaf", y = 2, label = "b", parse = TRUE, size=10) 

  lm <-  lm(pop0_mk_NA ~ log2FoldChangeAbs, interExp[interExp$sexBias == "Leaf",])
  summary(lm)  
  
  multiplot(X_mk,Y_mk,cols=2)
  
 interExp$pop0_mk_NA[interExp$log2FoldChangeAbsRound == 0 & interExp$sexBias == "Male"]

  
ggplot(interExp, aes(y=pop0_mk_NA, x=log2FoldChangeAbsRound, color=sexBias)) + 
  #geom_density_ridges() +
  #geom_violin() + 
  geom_boxplot(width=0.5) +
  geom_point(position=position_dodge(.5),aes(alpha=0.1) ) + guides(alpha=FALSE) +
  theme_bw()  + theme_bw(base_size = 30) + 
  #facet_grid(. ~ regions) +
  labs(x = "log Expression difference", y=title_mk,color="Bias") +
  scale_x_discrete(limits=c("0","1","2","3","4","5")) #+ ylim(0,5)

aov(pop0_mk_NA ~ log2FoldChangeAbsRound + sexBias + log2FoldChangeAbsRound*sexBias, data=interExp)


t.test(interExp$pop0_mk_NA[interExp$sexBias == "Unbiased"], interExp$pop0_mk_NA[interExp$sexBias == "Pollen"])
t.test(interExp$pop0_mk_NA[interExp$sexBias == "Leaf"],interExp$pop0_mk_NA[interExp$sexBias == "Pollen"])
t.test(interExp$pop0_mk_NA[interExp$sexBias == "Unbiased"],interExp$pop0_mk_NA[interExp$sexBias == "Leaf"])


mean(!is.na(interExp$pop0_mk_NA[interExp$sexBias == "Leaf"]))
mean(!is.na(interExp$pop0_mk_NA[interExp$sexBias == "Pollen"]))

multiplot(A_dnds,A_mk,cols=2)

##Fst##

ggplot(interExp, aes(x = log2FoldChange, y = pop0_Fst_syn, color=pvalue)) + 
  geom_point() + 
  stat_smooth(method = "loess") +
  # facet_grid(. ~ regions) +
  labs(x = "log2 Male Expression Bias", y="M-F Fst") +
  theme_bw()  + theme_bw(base_size = 30) + 
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    # '#1B75BB', #purple-y
    '#8B4BD8'
  ))

####ds####

interExp <- statsExp(set="rna",chrom="X2Y",highbias="Leaf",lowbias="Pollen",pop="phase")

ggplot(interExp, aes(y = log2FoldChange, x = pop0_d_syn, color=sexBias)) + 
geom_point() + 
  stat_smooth(method = "lm") +
  # facet_grid(. ~ regions) +
  labs(y = "log Leaf Expression Bias", x="X-Y ds",color="Tissue") +
  theme_bw()  + theme_bw(base_size = 30) 


lm <-  lm(pop0_d_syn ~ log2FoldChange, interExp[interExp$sexBias == "Pollen",])
summary(lm)  

head(allExp)





  
#sql bind FPKM with summarystats
  #instead of asking if genes are getting more or less biased, ask if they are more or less expressed


####relatedness between pollen samples####


####to genome####
interExp <- interExp[!is.na(interExp$pop0_pi_syn),]

interSplit <- separate(interExp, locus, c("locus","transcript","tail"), 
                       sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")

genomeConversion <- fread("NCF1_Genome.list",header=FALSE)

#genomeLength <- fread("Genome_length.txt",sep=" ")
#genome <- data.frame(sqldf('select genomeConversion.*, genomeLength.* from genomeConversion left join genomeLength on genomeConversion.V2 = genomeLength.V1'))
#genome <- genome[!is.na(genome$V1),]


#genome <- genome[genome$V2..14 > 500000,]
#sort(genome$V2..14)

S11619 <- genomeExp[genomeExp$Scaffold == "ScnbKXS_11619",]

genomeS <- separate(genomeConversion, V1, c("locus","transcript","tail"), 
                       sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")

genomeSplit <- separate(genomeS, V2, c("Scaffold","ScaffoldExtra"), 
                        sep = ";", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")


genomeExp <- data.frame(sqldf('select genomeSplit.*, interSplit.* from genomeSplit left join interSplit on interSplit.transcript = genomeSplit.transcript'))

genomeExp <- genomeExp[!is.na(genomeExp$pop0_pi_syn),]

ggplot(data=genomeExp) + 
  geom_point(aes(x=V9,y=pop1_pi_syn,color="Texas")) +
  geom_point(aes(x=V9,y=pop2_pi_syn,color="NC")) +
  facet_grid(aes(cols=Scaffold)) + 
  theme_bw()  + theme_bw(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  labs(x = "Scaffold Posittion", y=title_pisyn,color="Population") +
  scale_color_manual(values=c( 
    '#00ADEF', #Blue
    #'#FFF100',  #yellow
    #'#00A550', #green
    # '#1B75BB', #purple-y
    '#8B4BD8'
  ))
