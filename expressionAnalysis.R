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

base <- fread('allLoci.list',header = FALSE)
inds <- fread('rna.sql.list',header = FALSE)
inds <- inds$V1
inds <- c(inds,"SC1pollen","SC2pollen","TX1Bpollen","TX2Bpollen")

####MERGE####
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

readsCounts <- reads[,-c(1)]
FPKM <-  data.frame(
                    t(
                      t(readsCounts/
                          (readsCounts$length/1000)
                        )/(colSums(readsCounts, na.rm = FALSE, dims = 1)/1000000)
                      )
                    )


####DESEQ####
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

treads <- reads[,-c(1:2)]
rownames(treads) <- reads$V1

pop <- c("T","T","T","T","T","T","T","T","T","T","T","T","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","T","N","N","N","N","N","T","T","N","N","N","T","T","T","T","T","N","N","T","T")
sex <- c("M","M","M","F","F","F","F","F","F","M","M","M","F","F","F","F","F","F","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M")
origin <- c("J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","G","G","G","G")
tissue <- c("L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","P","P","P","P")
colData <- cbind(sex,pop,origin,tissue)
rownames(colData) <- colnames(reads[-c(1,2)])

countdata <- data.matrix(treads,rownames.force=TRUE) 

dds <- DESeqDataSetFromMatrix(countData = treads, colData = colData, design = ~  tissue )
dds <- dds[ rowSums(counts(dds)) > 1, ] #filter for rows with info

dss <- DESeq(dds)
#res <- results(dss, contrast=c("sex","M","F"),  alpha = 0.05  )
#resLFC <- lfcShrink(dss, coef="sex_M_vs_F", type="apeglm")
res <- results(dss, contrast=c("tissue","L","P"),  alpha = 0.05  )
resLFC <- lfcShrink(dss, coef="tissue_L_vs_P", type="apeglm")


head(res)
summary(res)
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))


####Orthologs####

roth <- fread('rothschildianus.orthologs',header=FALSE)
buc <- fread('bucephalophorus.orthologs',header=FALSE)

species <- data.frame(sqldf('select base.*, roth.*, buc.* from base left join roth on base.V1 = roth.V1 left join buc on base.V1 = buc.V1'))
species <- species[,c(1,3,5)]
names(species) <- c("hast","roth","buc")


var(between)/var(within) -> (tot - within) / tot

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
names(regions) <- c("C","A","XY","H","N")

regions$A[ !is.na(regions$A) == TRUE ] <- "A"
regions$XY[ !is.na(regions$XY) == TRUE ] <- "XY"
regions$H[ !is.na(regions$H) == TRUE ] <- "H"
regions$N[ !is.na(regions$N) == TRUE ] <- "N"

regions$regions[ !is.na(regions$A) == TRUE ] <- "Autosomal"
regions$regions[ !is.na(regions$XY) == TRUE ] <- "XY"
regions$regions[ !is.na(regions$H) == TRUE ] <- "Hemizygous"
regions$regions[ !is.na(regions$N) == TRUE ] <- "NeoXY"


resdf <- data.frame(res)
resdf$loci <- row.names(resdf)

regionsExp <- data.frame(sqldf('select resdf.*, regions.* from resdf left join regions on resdf.loci = regions.C'))

ggplot(regionsExp, aes(x = log2FoldChange, y = regions, fill=regions, alpha=0.5)) + 
  geom_density_ridges() + guides(alpha=FALSE,fill=FALSE) +
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "log2 Leaf Expression Bias", y="") +
scale_fill_manual(values=c( 
  '#00ADEF', #Blue
  '#FFF100',  #yellow
  '#00A550', #green
 # '#1B75BB', #purple-y
  '#8B4BD8'
))

#relatedness between pollen samples


inter <- fread('rna_rothschildianus_interpop_pop_A.txt')
within <- fread('rna_rothschildianus_summarystats_pop_A.txt')


inter <- separate(inter, locus, c("locus","file"), 
                    sep = ".fasta", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")

within <- separate(within, locus, c("locus","file"), 
                  sep = ".fasta", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")


interExp <- data.frame(sqldf('select regionsExp.*, inter.*, within.* from regionsExp left join inter on regionsExp.loci = inter.locus left join within on regionsExp.loci = within.locus'))

interExpA <- interExp[interExp$regions == "Autosomal" ,]

#interExpA$sexBias[interExpA$log2FoldChange > 0] <- "Male"
#interExpA$sexBias[interExpA$log2FoldChange < 0] <- "Female"
interExpA$sexBias[interExpA$log2FoldChange > 0] <- "Leaf"
interExpA$sexBias[interExpA$log2FoldChange < 0] <- "Pollen"
interExpA$sexBias[interExpA$pvalue >  0.05] <- "Unbiased"

pollenList <- interExpA$loci[interExpA$sexBias == "Pollen"]

write(pollenList, file = "pollen.list",
      append = FALSE, sep = "\n")

ggplot(interExpA, aes(x = log2FoldChange, y = pop0_tajD_syn, color=pvalue)) + 
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

ggplot(interExpA, aes(y=pop0_tajD_syn, x=sexBias)) + 
  geom_violin() +
  scale_x_discrete(limits=c("Leaf","Unbiased","Pollen")) +
  labs(x = "Expression Bias", y="TajD_syn") +
  theme_bw()  + theme_bw(base_size = 30)

t.test(interExpA$pop0_tajD_syn[interExpA$sexBias == "Leaf"],interExpA$pop0_tajD_syn[interExpA$sexBias == "Pollen"])

TajD_LP_anova <- aov(pop0_tajD_syn ~ sexBias, data=interExpA)
summary(TajD_LP_anova) 

