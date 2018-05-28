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

####ExtraFunctions####

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

summaryStatsInter <- function(inter=NULL,chrom=NULL){
  inter_melt <- melt(inter,id.vars = "locus")
  inter_comp <- inter_melt[complete.cases(inter_melt), ]
  inter_summary <- summarySE(inter_comp, measurevar="value", groupvars=c("variable"))
  cbind(inter_summary, "chrom" = chrom)
}

####import####

pops <- c("All","TX","NC")
wIn_roth_y <-
  fread('XYphased_rothschildianus_summarystats_pop_Ychrom.txt')
btw_roth_y <-
  fread('XYphased_rothschildianus_interpop_pop_Ychrom.txt')

XYph_r_Y_within <- summaryStats(in_mat=within,pops=pops,chrom="Y")
XYph_r_Y_inter <- summaryStatsInter(inter=between,chrom="Y")

####DXY####
XYph_r_Y_k <- XYph_r_Y_within[XYph_r_Y_within$var == "dxy" & XYph_r_Y_within$cod == "syn",]
XYph_r_Y_dxy <- XYph_r_Y_inter[XYph_r_Y_inter$variable == "Dxy_1_2",]

all_dxy <- rbind( XYph_r_Y_k[,-c(1,3)],
                  rename( XYph_r_Y_dxy,c("variable"="pop"))
                  )
all_dxy_out <- cbind(all_dxy,outgroup="rothschildianus")

ggplot(all_dxy_out, aes(x=outgroup, y=value, fill=pop)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Dxy") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  facet_grid(. ~ chrom, scales = "free")



