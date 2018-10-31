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
options(max.print = 1000000)

####import####

base <- fread('allLoci.list',header = FALSE)
inds <- fread('rna.sql.list',header = FALSE)
inds <- inds$V1
inds <- c(inds,"SC1pollen","SC2pollen","TX1Bpollen","TX2Bpollen")

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

#readsCounts <- reads[,-c(1)]
#FPKM <-  data.frame(t( t(readsCounts/  (readsCounts$length/1000)  )/(colSums(readsCounts, na.rm = FALSE, dims = 1)/1000000) ) )

####DESEQ####
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#Sample Info#
pop <- c("T","T","T","T","T","T","T","T","T","T","T","T","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","T","N","N","N","N","N","T","T","N","N","N","T","T","T","T","T","N","N","T","T")
sex <- c("M","M","M","F","F","F","F","F","F","M","M","M","F","F","F","F","F","F","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M")
origin <- c("J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","G","G","G","G")
tissue <- c("L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","P","P","P","P")
colData <- cbind(sex,pop,origin,tissue)
rownames(colData) <- colnames(reads[-c(1,2)])

#Read Intake#
treads <- reads[,-c(1:2)]
rownames(treads) <- reads$V1
countdata <- data.matrix(treads,rownames.force=TRUE) 

##Define bias analysis##
#dds <- DESeqDataSetFromMatrix(countData = treads, colData = colData, design = ~ tissue)
dds <- DESeqDataSetFromMatrix(countData = treads, colData = colData, design = ~  sex + tissue)

dds <- dds[ rowSums(counts(dds)) > 1, ] #filter for rows with info
dss <- DESeq(dds)

res <- results(dss, contrast=c("sex","M","F"),  alpha = 0.05  )
#resLFC <- lfcShrink(dss, coef="sex_M_vs_F", type="apeglm")

#res <- results(dss, contrast=c("tissue","L","P"),  alpha = 0.05  )
#resLFC <- lfcShrink(dss, coef="tissue_L_vs_P", type="apeglm")

#plotMA(res, ylim=c(-2,2))
#plotMA(resLFC, ylim=c(-2,2))

resdf <- data.frame(res)
resdf$loci <- row.names(resdf)

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
regionsExp <- regionsExp[regionsExp$pvalue <  0.05,]
regionsExp <- regionsExp[!is.na(regionsExp$pvalue), ]


ggplot(regionsExp, aes(x = log2FoldChange, y = regions, fill=regions, alpha=0.5)) + 
  geom_density_ridges() + guides(alpha=FALSE,fill=FALSE) +
  theme_bw()  + theme_bw(base_size = 30) + 
  labs(x = "log2 Male Expression Bias", y="") +
 # labs(x = "log2 Leaf Expression Bias", y="") +
scale_fill_manual(values=c( 
  '#00ADEF', #Blue
  '#FFF100',  #yellow
  '#00A550', #green
 # '#1B75BB', #purple-y
  '#8B4BD8'
)) + 
  xlim(-3,3)




####Stat Correlates####

inter <- fread('rna_rothschildianus_interpop_fm_A.txt')
within <- fread('rna_rothschildianus_summarystats_fm_A.txt')

inter <- separate(inter, locus, c("locus","file"), 
                    sep = ".fasta", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")

within <- separate(within, locus, c("locus","file"), 
                  sep = ".fasta", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")

interExp <- data.frame(sqldf('select regionsExp.*, inter.*, within.* from regionsExp 
                             left join inter on regionsExp.loci = inter.locus 
                             left join within on regionsExp.loci = within.locus'))

#interExp <- interExp[interExp$regions == "Autosomal" ,]

interExp$sexBias[interExp$log2FoldChange > 0] <- "Male"
interExp$sexBias[interExp$log2FoldChange < 0] <- "Female"
#interExp$sexBias[interExp$log2FoldChange > 0] <- "Leaf"
#interExp$sexBias[interExp$log2FoldChange < 0] <- "Pollen"
interExp$sexBias[interExp$pvalue >  0.05] <- "Unbiased"

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

stacks <- data.frame(
cbind(
  rbind(
    length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]) / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ])),
    length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]) / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'XY' ]))
  ),
  rbind(
    length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ] ) / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ])),
    length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ] ) / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Autosomal' ]))
  ),
  rbind(
    length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ] ) / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ])),
    length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ] ) / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'Hemizygous' ]))
  ),
  rbind(
    length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ] ) / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ])),
    length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ] ) / (  length(interExp$log2FoldChange[interExp$log2FoldChange < 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ]) +length(interExp$log2FoldChange[interExp$log2FoldChange > 0 & interExp$pvalue <  0.05 & interExp$regions == 'NeoXY' ]))
  )
)
,row.names = c("Female","Male"))
#,row.names = c("Pollen","Leaf"))


names(stacks) <- c("XY","Autosomal","Hemizygous","NeoXY")
tstacks <- data.frame(t(stacks))
tstacks$region <- rownames(tstacks)
stack_melt <- melt(tstacks,id.vars = "region")

ggplot(stack_melt, aes(x = region, y = value,fill=variable)) + 
  geom_bar(stat="identity") +
  # facet_grid(. ~ regions) +
  labs(x = "", y="Percent Expression Bias",fill="") +
  theme_bw()  + theme_bw(base_size = 30) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    #'#FFF100',  #yellow
    #'#00A550', #green
    # '#1B75BB', #purple-y
    '#8B4BD8'
  )) 



#Tajima's D#

ggplot(interExp, aes(x = log2FoldChange, y = pop0_tajD_syn, color=pvalue)) + 
  geom_point() + 
  stat_smooth(method = "loess") +
 # facet_grid(. ~ regions) +
  labs(x = "log2 Male Expression Bias", y="TajD_syn") +
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
  labs(x = "Expression Bias", y=title_tds) +
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
  labs(x = "Expression Bias", y=title_pisyn) +
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



##

ggplot(interExpBiased_comp, aes(y=pop0_dnds_NA, x=log2FoldChangeAbs, color=sexBias)) + 
  geom_point() + guides(alpha=FALSE) +
  stat_smooth(method = "lm",se=FALSE,size=2) +
  theme_bw()  + theme_bw(base_size = 30) +
  labs(x = "Expression Bias", y="dnds", color="Sex") 

#Fst

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

####Orthologs####

roth <- fread('rothschildianus.orthologs',header=FALSE)
buc <- fread('bucephalophorus.orthologs',header=FALSE)

species <- data.frame(sqldf('select base.*, roth.*, buc.* from base left join roth on base.V1 = roth.V1 left join buc on base.V1 = buc.V1'))
species <- species[,c(1,3,5)]
names(species) <- c("hast","roth","buc")


var(between)/var(within) -> (tot - within) / tot

####relatedness between pollen samples####