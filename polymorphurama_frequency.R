
setwd('~/Google Drive/Research/Data/')
library(data.table)
library(ggplot2)
#import .frequency

within_stat <- fread('fem_bucephalophorus_summarystats_pop_H.txt' )

invar_syn <- sum(within_stat$pop1_sites_syn  - within_stat$pop1_S_syn )
invar_rep <- sum( within_stat$pop1_sites_rep - within_stat$pop1_S_rep)


data <- fread('fem_bucephalophorus_frequencies_pop_H.txt',header = FALSE )
length <- ncol(data) / 2
cat(length)

data_syn <- data[,2:length]
data_rep <- data[,(length+2):(length*2)]

data_syn_tot <- colSums (data_syn, na.rm = FALSE, dims = 1)
cat(data_syn_tot)

data_rep_tot <- colSums (data_rep, na.rm = FALSE, dims = 1)
cat(invar_rep,data_rep_tot)


data_syn_rel<-data.frame((data_syn_tot/sum(data_syn_tot)))

ggplot(data_syn_rel) + geom_histogram(aes(x=X.data_syn_tot.sum.data_syn_tot..)) +
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "Count", y="No. alleles") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) 
  
data_rep_rel<-data.frame((data_rep_tot/sum(data_rep_tot)))

ggplot(data_rep_rel) + geom_histogram(aes(x=X.data_rep_tot.sum.data_rep_tot..)) +
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "Count", y="No. alleles") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) 

ggplot() + 
  geom_histogram(aes(x=data_syn_rel$X.data_syn_tot.sum.data_syn_tot..,alpha=0.1,fill="syn")) +
  geom_histogram(aes(x=data_rep_rel$X.data_rep_tot.sum.data_rep_tot..,alpha=0.1,fill="rep")) +
  theme_bw()  + theme_bw(base_size = 30) + labs(y = "Count", x="No. alleles") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
  guides(alpha=FALSE,fill=guide_legend(title="Sites")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100'
  ))



