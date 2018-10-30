
#library(ggplot2)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
#setwd('~/Google Drive/Research/Data/')

set= args[1]
outgroup=args[2]
popStr=args[3]
popNum=args[4]
chrom=args[5]

filename_wIn <- paste(set,"_",outgroup,"_summarystats_",popStr,"_",chrom,".txt",sep="")
#filename_wIn <- "rna_rothschildianus_summarystats_pop_pollen.txt"
within_stat <- fread(filename_wIn)
invar_syn <- sum(within_stat$pop1_sites_syn  - within_stat$pop1_S_syn )
invar_rep <- sum( within_stat$pop1_sites_rep - within_stat$pop1_S_rep)

#filename_frq <- "rna_rothschildianus_frequencies_pop1_pollen.txt"
filename_frq <- paste(set,"_",outgroup,"_frequencies_",popStr,"1_",chrom,".txt",sep="")
data <- fread(filename_frq,header = FALSE )
length <- (ncol(data) / 2)-1
cat("1\n")
cat(length,"\n")

data_syn <- data[,2:(length+1)]
data_rep <- data[,(length+3):(length*2+2)]

data_syn_tot <- colSums (data_syn, na.rm = FALSE, dims = 1)
data_rep_tot <- colSums (data_rep, na.rm = FALSE, dims = 1)

data_rep_fold <- vector("list", (length+1))
data_rep_fold[1] <- round(invar_rep)

data_rep_fold[2] <- data_rep_tot[2]
for (site in c(seq(1,(length/2),by=1))){
  data_rep_fold[site+2] <-  data_rep_tot[site+2] + data_rep_tot[(length+1)-site]
  data_rep_fold[length+2-site] <- 0
}


#for (site in c(seq(1,ceiling(length/2),by=1))){
#  data_rep_fold[site+1] <-  data_rep_tot[site] + data_rep_tot[(length+1)-site]
#  data_rep_fold[length+2-site] <- 0
#}


do.call(cat,data_rep_fold)
cat("\n")

data_syn_fold <- vector("list", (length+1))
data_syn_fold[1] <- round(invar_syn)
data_syn_fold[2] <- data_syn_tot[2]

for (site in c(seq(1,((length)/2),by=1))){
  data_syn_fold[site+2] <-  data_syn_tot[site+2] + data_syn_tot[(length+1)-site]
  data_syn_fold[length+2-site] <- 0
}

#for (site in c(seq(1,ceiling(length/2),by=1))){
#  data_syn_fold[site+1] <-  data_syn_tot[site] + data_syn_tot[(length+1)-site]
#  data_syn_fold[length+2-site] <- 0
#}

do.call(cat,data_syn_fold)
cat("\n")
  
#data_syn_rel<-data.frame((data_syn_tot/sum(data_syn_tot)))

#ggplot(data_syn_rel) + geom_histogram(aes(x=X.data_syn_tot.sum.data_syn_tot..)) +
#  theme_bw()  + theme_bw(base_size = 30) + labs(x = "Count", y="No. alleles") +
#  theme(axis.text.x = element_text(angle = 40, hjust = 1)) 
  
#data_rep_rel<-data.frame((data_rep_tot/sum(data_rep_tot)))

#ggplot(data_rep_rel) + geom_histogram(aes(x=X.data_rep_tot.sum.data_rep_tot..)) +
#  theme_bw()  + theme_bw(base_size = 30) + labs(x = "Count", y="No. alleles") +
#  theme(axis.text.x = element_text(angle = 40, hjust = 1)) 

#ggplot() + 
#  geom_histogram(aes(x=data_syn_rel$X.data_syn_tot.sum.data_syn_tot..,alpha=0.1,fill="syn")) +
#  geom_histogram(aes(x=data_rep_rel$X.data_rep_tot.sum.data_rep_tot..,alpha=0.1,fill="rep")) +
#  theme_bw()  + theme_bw(base_size = 30) + labs(y = "Count", x="No. alleles") +
#  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
#  guides(alpha=FALSE,fill=guide_legend(title="Sites")) +
#  scale_fill_manual(values=c( 
#    '#00ADEF', #Blue
 #   '#FFF100'
#  ))


#SFS <- fread('pollen_SFS_OE.txt')

#SFS_melt <- melt(SFS,id.vars = "Alleles",verbose=FALSE)

#ggplot(SFS_melt,aes(x=Alleles, y=value,fill=variable)) +
#  geom_bar(stat="identity", color="black", position=position_dodge())+
#  theme_minimal() + labs(x = "Count", y="Frequency",fill="")

#dfe <- fread('dfe_pollen.txt')
#dfe$cat[dfe$V1 == 0] <- "0<Nes<1"
#dfe$cat[dfe$V1 == 1] <- "1<Nes<10"
#dfe$cat[dfe$V1 == 10] <- "10<Nes<100"
#dfe$cat[dfe$V1 == 100] <- "100<Nes"


#ggplot(dfe,aes(x=cat, y=V2)) +
#  geom_bar(stat="identity", color="black", position=position_dodge())+
#  theme_minimal() +  theme_bw(base_size = 30) +
#  labs(x = "", y="Proportion",title="Texas Pollen-Bias") 

