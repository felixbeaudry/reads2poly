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
    pol_melt <- melt(in_mat,id.vars = "locus",verbose=FALSE)
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
    inter_melt <- melt(inter,id.vars = "locus",verbose=FALSE)
    inter_sep <- separate(inter_melt, variable, c("pop","var","cod"), 
                          sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    inter_comp <- inter_sep[complete.cases(inter_sep), ]
    
    for (i in 1:(length(pops))){
      popCount = paste("pop",(i-1),sep="") 
      inter_comp$pop[inter_comp$pop == popCount] <- pops[i]
    }
    
    inter_summary <- summarySE(inter_comp, measurevar="value", groupvars=c("var","pop","cod"))
    cbind(inter_summary, "chrom" = chrom)
  }
  filename_wIn <- paste(set,"_",outgroup,"_summarystats_pop_",chrom,"chrom.txt",sep="")
  filename_btw <- paste(set,"_",outgroup,"_interpop_pop_",chrom,"chrom.txt",sep="")

  in_read <- fread(filename_wIn)
  seqMax <- max(in_read$pop0_seqs_NA[!is.na(in_read$pop0_seqs_NA)])
  coverage_list <- in_read$locus[in_read$pop0_seqs_NA == seqMax]
  in_read_cov <- in_read[in_read$locus %in% coverage_list ]
  in_bet <- fread(filename_btw)
  in_bet_cov <- in_bet[in_bet$locus %in% coverage_list ]
  
  stats_wIn <- summaryStats(in_mat=in_read_cov ,pops=pops,chrom=chrom)
  stats_btw <- summaryStatsInter(inter=in_bet_cov,pops=pops,chrom=chrom)
  stats_bind <- cbind(rbind(stats_wIn,stats_btw),outgroup=outgroup,stringsAsFactors=FALSE)
  return(stats_bind)
}
stats_var <- function(outgroup=NULL,set=NULL,chrom=NULL,pops=NULL){ 

  summaryVar <- function(in_mat=NULL,pops=NULL,chrom=NULL){
    pol_melt <- melt(in_mat,id.vars = "locus")
    pol_sep <- separate(pol_melt, variable, c("pop","var","cod"), 
                        sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    pol_sep_comp <- pol_sep[complete.cases(pol_sep), ]
    
    for (i in 1:(length(pops))){
      popCount = paste("pop",(i-1),sep="") 
      pol_sep_comp$pop[pol_sep_comp$pop == popCount] <- pops[i]
    }
    cbind(pol_sep_comp, "chrom" = chrom)
  }
  
  
  summaryVarInter <- function(inter=NULL,pops=NULL,chrom=NULL){
    inter_melt <- melt(inter,id.vars = "locus")
    inter_sep <- separate(inter_melt, variable, c("pop","var","cod"), 
                          sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    inter_comp <- inter_sep[complete.cases(inter_sep), ]
    
    for (i in 1:(length(pops))){
      popCount = paste("pop",(i-1),sep="") 
      inter_comp$pop[inter_comp$pop == popCount] <- pops[i]
    }
    
    cbind(inter_comp, "chrom" = chrom)
  }
  filename_wIn <- paste(set,"_",outgroup,"_summarystats_pop_",chrom,"chrom.txt",sep="")
  filename_btw <- paste(set,"_",outgroup,"_interpop_pop_",chrom,"chrom.txt",sep="")
  stats_wIn <- summaryVar(in_mat=fread(filename_wIn),pops=pops,chrom=chrom)
  stats_btw <- summaryVarInter(inter=fread(filename_btw),pops=pops,chrom=chrom)
  stats_bind <- cbind(rbind(stats_wIn,stats_btw),outgroup=outgroup,stringsAsFactors=FALSE)
  return(stats_bind)
}
ms_stat <- function(chrom=NULL,var=NULL){
  filename <- paste('ms_',chrom,'_stat.txt',sep="")
  ms <- fread(filename)
  ms_melt <- melt(ms,id.vars = "rep")
  ms_comp <- ms_melt[complete.cases(ms_melt), ]
  ms_sum <- summarySE(ms_comp, measurevar="value", groupvars="variable")
  #chrom_name <- paste(chrom,'_sim',sep="")
  ms_bind <- cbind(rename(ms_sum, c("variable" = "var")),"chrom"=chrom,"outgroup"="simulation","state"="simulation")
  ms_out <- ms_bind[ms_bind$var == var,]
  return(ms_out)
}


####import####

#outgroup <- c("rothschildianus","bucephalophorus")
#set <- c("rna","rna_felix","rna_josh","rna_josh_males","rna_josh_fem")
#chrom <- c("auto","hemi","X","Y")

#pop <- c("hastatulus","Y1","Y1Y2")
pop <- c("Rhastatulus","XY","XYY")

all_data <- data.frame(rbind(
stats_table(outgroup="rothschildianus",set="XYphased",chrom="X",pops=pop)
,stats_table(outgroup="rothschildianus",set="XYphased",chrom="Y",pops=pop)
,stats_table(outgroup="rothschildianus",set="RNA",chrom="H",pops=pop)
,stats_table(outgroup="bucephalophorus",set="RNA",chrom="H",pops=pop)
,stats_table(outgroup="bucephalophorus",set="RNA",chrom="A",pops=pop)
,stats_table(outgroup="bucephalophorus",set="XYphased",chrom="X",pops=pop)
,stats_table(outgroup="bucephalophorus",set="XYphased",chrom="Y",pops=pop)
), stringsAsFactors = FALSE)

all_data_var<- data.frame(rbind(
  stats_var(outgroup="rothschildianus",set="XYphased",chrom="X",pops=pop)
  ,stats_var(outgroup="rothschildianus",set="XYphased",chrom="Y",pops=pop)
  ,stats_var(outgroup="rothschildianus",set="RNA",chrom="H",pops=pop)
  ,stats_var(outgroup="bucephalophorus",set="RNA",chrom="H",pops=pop)
  ,stats_var(outgroup="bucephalophorus",set="RNA",chrom="A",pops=pop)
  ,stats_var(outgroup="bucephalophorus",set="XYphased",chrom="X",pops=pop)
  ,stats_var(outgroup="bucephalophorus",set="XYphased",chrom="Y",pops=pop)
), stringsAsFactors = FALSE)

maleX <- data.frame(stats_table(outgroup="rothschildianus",set="XYphased",chrom="X",pops=pop))
femaleX <- data.frame(stats_table(outgroup="rothschildianus",set="joshFem",chrom="X",pops=pop))

maleX[maleX$var == "pi",]
femaleX[maleX$var == "pi",]

####pi####

all_data_pi <- all_data[all_data$var == "pi" & all_data$cod == "syn"  ,]

title_pisyn <- expression(paste(pi, ""[syn]))

ggplot(all_data_pi, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_pisyn) +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  facet_grid(. ~ pop, scales = "free") +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#1B75BB', #purple-y
    '#00A550' #green
  ))


ms_pi <- rbind(ms_stat(chrom="X",var="pi_tot"),
                ms_stat(chrom="A",var="pi_tot"),
               all_data_pi[all_data_pi$pop == "Rhastatulus",-c(2,3)]
)

ggplot(ms_pi, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_pisyn) +
  scale_x_discrete(limits=c("A","A_sim","H","X","X_sim","Y")) 

####TajD####

all_data_tajD <- all_data[all_data$var == "tajD" & all_data$cod == "syn"  ,]

ggplot(all_data_tajD, aes(x=chrom, y=value, fill=chrom)) + 
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Tajima's D") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  facet_grid(. ~ pop, scales = "free") +
  
  scale_x_discrete(limits=c("A","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#1B75BB', #purple-y
    '#00A550' #green
  ))


all_data[all_data$var == "sites" && all_data$pop == "Rhastatulus",]

####FST####

all_data_fst <- all_data[all_data$var == "Fst" & all_data$cod == "syn"  ,]



ggplot(all_data_fst, aes(x=chrom, y=value, fill=chrom)) +
 guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Fst") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#1B75BB', #purple-y
    '#00A550' #green
  ))

ms_fst <- rbind(ms_stat(chrom="X",var="fst"),
                ms_stat(chrom="A",var="fst"),
                cbind(all_data_fst[,-c(2,3)],"state"="empirical")
)

#fst_plot <- 
  ggplot(ms_fst, aes(x=chrom, y=value, fill=state)) +
 guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Fst") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  scale_x_discrete(limits=c("A","H","X","Y")) 

####dxy####

dxy <- all_data[all_data$var == "dxy" ,]

ggplot(dxy, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Dxy") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#1B75BB', #purple-y
    '#00A550' #green
    ))


ms_dxy <- rbind(
  ms_stat(chrom="X",var="dxy"),
  ms_stat(chrom="A",var="dxy"),
  cbind(dxy[,-c(2,3)],"state"="empirical")
)

#dxy_plot <- 
  ggplot(ms_dxy, aes(x=chrom, y=value, fill=state)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Dxy") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  scale_x_discrete(limits=c("A","H","X","Y")) 

#multiplot(fst_plot,dxy_plot,cols=2)

t.test(
  all_data_var$value[all_data_var$chrom == "A" & all_data_var$var == "dxy" ] ,  
  all_data_var$value[all_data_var$chrom == "Y" & all_data_var$var == "dxy" ])

dxy_data_var <- all_data_var[ all_data_var$var == "dxy", ]  
dxy_anova <- aov(value ~ chrom, data=dxy_data_var)
summary(dxy_anova) 


####kxy####

all_data_dxy <- all_data[all_data$var == "kxy" | all_data$var == "dxy",]
#levels(all_data_dxy$outgroup) = c("rothschildianus", "hastatulus","bucephalophorus")
all_data_dxy$outgroup[all_data_dxy$var == "dxy"] <- "hastatulus"


ggplot(all_data_dxy, aes(x=outgroup, y=value, color=chrom
                      # , group=1
                       )) + #guides(fill = FALSE) +
  geom_point(position=position_dodge(.9)) +
  
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Species",y="Average Divergence") +
  scale_x_discrete(limits=c("hastatulus","rothschildianus","bucephalophorus"))

####dnds####

evo_rate <- all_data[all_data$var == "dnds" | all_data$var == "kaks",]
#levels(all_data_dxy$outgroup) = c("rothschildianus", "hastatulus","bucephalophorus")
evo_rate$outgroup[evo_rate$var == "dnds"] <- "hastatulus"


ggplot(evo_rate, aes(x=outgroup, y=value, color=chrom
                         # , group=1
)) + #guides(fill = FALSE) +
  geom_point(position=position_dodge(.9)) +
  
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Species",y="Average Rate of Evolution") +
  scale_x_discrete(limits=c("hastatulus","rothschildianus","bucephalophorus"))

####alpha####

alpha <- all_data[all_data$var == "alpha" ,]

alpha$chrom_f = factor(alpha$chrom, levels=c('A','H','X','Y'))

title_alpha <- expression(paste(alpha))

ggplot(alpha, aes(x=pop, y=value, color=outgroup
                     # , group=1
)) + #guides(fill = FALSE) +
  geom_point(position=position_dodge(.9)) +
  
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="Population",y=title_alpha) +
  facet_grid(. ~ chrom_f, scales = "free") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))  +
  theme(strip.background =element_rect(fill="white"))
