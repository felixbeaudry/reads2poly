#### TOP ####
setwd('~/Google Drive/Research/Data/')
library(data.table)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(sqldf)
options(max.print = 1000000)

####HomemadeFunctions####

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

stats_table <- function(outgroup=NULL,set=NULL,chrom=NULL,pops=NULL){ 
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    library(plyr)
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )
    datac <- rename(datac, c("mean" = measurevar))
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    return(datac)
  }
  
  summaryStats <- function(in_mat=NULL,pops=NULL,chrom=NULL){
    pol_melt <- melt(in_mat,id.vars = "locus")
    pol_sep <- separate(pol_melt, variable, c("pop","var","cod"), 
                        sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    pol_sep_comp <- pol_sep[complete.cases(pol_sep), ]
    
    for (i in 1:(length(pops))){
      popCount = paste("pop",(i-1),sep="") 
      pol_sep_comp$pop[pol_sep_comp$pop == popCount] <- pops[i]
    }
    summary <- summarySE(pol_sep_comp, measurevar="value", groupvars=c("var","pop","cod"))
    cbind(summary, "chrom" = chrom)
  }
  
  
  summaryStatsInter <- function(inter=NULL,pops=NULL,chrom=NULL){
    inter_melt <- melt(inter,id.vars = "locus")
    inter_sep <- separate(inter_melt, variable, c("pop","var","cod"), 
                          sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    inter_comp <- inter_sep[complete.cases(inter_sep), ]
    
    for (i in 1:(length(pops))){
      popCount = paste("pop",(i-1),sep="") 
      pol_sep_comp$pop[pol_sep_comp$pop == popCount] <- pops[i]
    }
    
    inter_summary <- summarySE(inter_comp, measurevar="value", groupvars=c("var","pop","cod"))
    cbind(inter_summary, "chrom" = chrom)
  }
  filename_wIn <- paste(set,"_",outgroup,"_summarystats_pop_",chrom,"chrom.txt",sep="")
  filename_btw <- paste(set,"_",outgroup,"_interpop_pop_",chrom,"chrom.txt",sep="")
  stats_wIn <- summaryStats(in_mat=fread(filename_wIn),pops=pops,chrom=chrom)
  stats_btw <- summaryStatsInter(inter=fread(filename_btw),pops=pops,chrom=chrom)
  stats_bind <- cbind(rbind(stats_wIn,stats_btw),outgroup=outgroup)
  return(stats_bind)
}

####import####

#outgroup <- c("rothschildianus","bucephalophorus")
#set <- c("rna","rna_felix","rna_josh","rna_josh_males","rna_josh_fem")
#chrom <- c("auto","hemi","X","Y")

pops <- c("hastatulus","Y1","Y1Y2")

xyphased <- rbind(
stats_table(outgroup="rothschildianus",set="XYphased",chrom="X",pops=pops),
stats_table(outgroup="rothschildianus",set="XYphased",chrom="Y",pops=pops),
stats_table(outgroup="bucephalophorus",set="XYphased",chrom="X",pops=pops),
stats_table(outgroup="bucephalophorus",set="XYphased",chrom="Y",pops=pops)
)


####charts####

ggplot(all_dxy, aes(x=outgroup, y=value, fill=pop)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Dxy") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  facet_grid(. ~ chrom, scales = "free")

ggplot(all_dxy, aes(x=outgroup, y=value, colour=chrom)) + #guides(fill = FALSE) +
  geom_point() +
  geom_line(aes()) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Dxy") +
  scale_x_discrete(limits=c("hastatulus","rothschildianus","bucephalophorus"))

