#### TOP ####
setwd('~/Google Drive/Research/Data/')
library(data.table)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(sqldf)
library(rapport)
library(Hmisc)

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

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {

  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,.fun = function(xx, col) {
                   c(
                     N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 }, measurevar)
  
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

stats_table <- function(outgroup=NULL,set=NULL,chrom=NULL,pops=NULL,subsetList=NULL,popStr='pop',ksyn=1){ 

  ##Calculate summary for within population statistics
  summaryStats <- function(in_mat=NULL,pops=NULL,chrom=NULL,set=NULL){
    pol_melt <- melt(in_mat,id.vars = "locus",verbose=FALSE)
    pol_sep <- separate(pol_melt, variable, c("pop","var","cod"), 
                        sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    pol_sep_comp <- pol_sep[complete.cases(pol_sep), ]
    
    for (i in 1:(length(pops))){
      popCount = paste("pop",(i-1),sep="") 
      pol_sep_comp$pop[pol_sep_comp$pop == popCount] <- pops[i]
    }
    summary <- summarySE(pol_sep_comp, measurevar="value", groupvars=c("var","pop","cod"))
    cbind(summary, "chrom" = chrom, "sex" = set)

  }
  
  ##Calculate summary for between population statistics
  summaryStatsInter <- function(inter=NULL,pops=NULL,chrom=NULL,set=NULL){
    inter_melt <- melt(inter,id.vars = "locus",verbose=FALSE)
    inter_sep <- separate(inter_melt, variable, c("pop","var","cod"), 
                          sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    inter_comp <- inter_sep[complete.cases(inter_sep), ]
    
    for (i in 1:(length(pops))){
      popCount = paste("pop",(i-1),sep="") 
      inter_comp$pop[inter_comp$pop == popCount] <- pops[i]
    }
    
    inter_summary <- summarySE(inter_comp, measurevar="value", groupvars=c("var","pop","cod"))
    cbind(inter_summary, "chrom" = chrom, "sex" = set)
  }

  ##import files
  filename_wIn <- paste(set,"_",outgroup,"_summarystats_",popStr,"_",chrom,".txt",sep="")
  filename_btw <- paste(set,"_",outgroup,"_interpop_",popStr,"_",chrom,".txt",sep="")
  in_read <- fread(filename_wIn)
  in_bet <- fread(filename_btw)
  
  ##filter for min 80% of inds and greater than 50 sites
  seqMax <- max(in_read$pop0_seqs_NA[!is.na(in_read$pop0_seqs_NA)])
  coverage_list <- in_read$locus[in_read$pop0_seqs_NA >= (seqMax * 0.8)  & in_read$pop0_sites_syn >= 50 & in_read$pop0_k_syn < ksyn]
  
  in_read_cov <- in_read[in_read$locus %in% coverage_list ]
  in_bet_cov <- in_bet[in_bet$locus %in% coverage_list ]
  
  ##allow locus subsetting inline##
  #if (length(subsetList) >0 ){
  #  in_read_cov <- in_read[in_read$locus %in% subsetList ]
  #  in_bet_cov <- in_bet_cov[in_bet_cov$locus %in% subsetList ]
  #}
  
  ##intake and bind
  stats_wIn <- summaryStats(in_mat=in_read_cov ,pops=pops,chrom=chrom,set=set)
  stats_btw <- summaryStatsInter(inter=in_bet_cov,pops=pops,chrom=chrom,set=set)
  stats_bind <- cbind(rbind(stats_wIn,stats_btw),outgroup=outgroup,stringsAsFactors=FALSE)
  return(stats_bind)
}
stats_var <- function(outgroup=NULL,set=NULL,chrom=NULL,pops=NULL,popStr='pop',ksyn=1){ 

  summaryVar <- function(in_mat=NULL,pops=NULL,chrom=NULL){
    pol_melt <- melt(in_mat,id.vars = "locus")
    pol_sep <- separate(pol_melt, variable, c("pop","var","cod"), 
                        sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    pol_sep_comp <- pol_sep[complete.cases(pol_sep), ]
    
    for (i in 1:(length(pops))){
      popCount = paste("pop",(i-1),sep="") 
      pol_sep_comp$pop[pol_sep_comp$pop == popCount] <- pops[i]
    }
    cbind(pol_sep_comp, "chrom" = chrom, "sex" = set)
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
    
    cbind(inter_comp, "chrom" = chrom, "sex" = set)
  }
  filename_wIn <- paste(set,"_",outgroup,"_summarystats_",popStr,"_",chrom,".txt",sep="")
  filename_btw <- paste(set,"_",outgroup,"_interpop_",popStr,"_",chrom,".txt",sep="")
  
  in_read <- fread(filename_wIn)
  seqMax <- max(in_read$pop0_seqs_NA[!is.na(in_read$pop0_seqs_NA)])
  coverage_list <- in_read$locus[in_read$pop0_seqs_NA >= (seqMax * 0.8) & in_read$pop0_sites_syn >= 50 & in_read$pop0_k_syn < ksyn]
  in_read_cov <- in_read[in_read$locus %in% coverage_list ]
  in_bet <- fread(filename_btw)
  in_bet_cov <- in_bet[in_bet$locus %in% coverage_list ]
  
  
  stats_wIn <- summaryVar(in_mat=fread(filename_wIn),pops=pops,chrom=chrom)
  stats_btw <- summaryVarInter(inter=fread(filename_btw),pops=pops,chrom=chrom)
  stats_bind <- cbind(rbind(stats_wIn,stats_btw),outgroup=outgroup,stringsAsFactors=FALSE)
  return(stats_bind)
}

ms_stat <- function(chrom=NULL,var=NULL,sitemean=NULL,filename=filname){
  #filename <- paste('ms_',chrom,'_stat.txt',sep="")
  ms <- fread(filename)
  ms_melt <- melt(ms,id.vars = "rep")
  ms_comp <- ms_melt[complete.cases(ms_melt), ]
  ms_sum <- summarySE(ms_comp, measurevar="value", groupvars="variable")
  #chrom_name <- paste(chrom,'_sim',sep="")
  ms_bind <- cbind(rename(ms_sum, c("variable" = "var")),"chrom"=chrom,"outgroup"="exp","state"="exp")
  ms_out <- ms_bind[ms_bind$var == var,]
  if (var == "pi_tot" | var == "dxy"){
    ms_out$value <- ms_out$value / sitemean
  }
  return(ms_out)
}

t.testloop <- function(chromSet=NULL,var=NULL,pop=NULL,cod=NULL,outgroup=NULL){
  for (i in 1:(length(chromSet))){
    for (j in 1:(length(chromSet))){
      
      result <- t.test(
        all_data_var$value[ all_data_var$chrom == chromSet[i] & all_data_var$var == var & all_data_var$pop == pop & all_data_var$cod == cod & all_data_var$outgroup == outgroup  ]
        ,
        all_data_var$value[ all_data_var$chrom == chromSet[j] & all_data_var$var == var  & all_data_var$pop == pop & all_data_var$cod == cod & all_data_var$outgroup == outgroup  ]
      )
      if(result$p.value < 0.05){
        print(c(chromSet[i],chromSet[j],"pvalue=",result$p.value,'*'))
      }

      else if(result$p.value != 1){
        print(c(chromSet[i],chromSet[j],"pvalue=",result$p.value))
      }
    }
  }
}
  

####import####

pop <- c("R.hastatulus","XY","XYY")
mpop = fpop = pop
FLNC <- c("FLNC","FL","NC")

TXNC <- c("TXNC","TX","NC")
TXFL <- c("TXFL","TX","FL")
phase = X2Y = c("X2Y","X2","2Y") 

 all_data <- data.frame(
  rbind(
    #X Male 
    stats_table(outgroup="bucephalophorus",set="sub4",chrom="X",pops=pop,popStr="pop",ksyn = 0.4),
    stats_table(outgroup="rothschildianus",set="sub4",chrom="X",pops=pop,popStr="pop",ksyn = 0.3),
    #Y Male
    stats_table(outgroup="bucephalophorus",set="sub4",chrom="Y",pops=pop,popStr="pop",ksyn = 0.4),
    stats_table(outgroup="rothschildianus",set="sub4",chrom="Y",pops=pop,popStr="pop",ksyn = 0.3),
    #H Male
    stats_table(outgroup="rothschildianus",set="m",chrom="H",pops=pop,popStr="pop"),
   stats_table(outgroup="bucephalophorus",set="m",chrom="H",pops=pop,popStr="pop"),
    #A Male
  #  stats_table(outgroup="bucephalophorus",set="m",chrom="A",pops=pop,popStr="pop"),
   stats_table(outgroup="rothschildianus",set="m",chrom="A",pops=pop,popStr="pop"),
  
  stats_table(outgroup="rothschildianus",set="m",chrom="Y",pops=FLNC,popStr="FLNC"),
  stats_table(outgroup="rothschildianus",set="m",chrom="A",pops=FLNC,popStr="FLNC"),
  stats_table(outgroup="rothschildianus",set="m",chrom="H",pops=FLNC,popStr="FLNC"),
  stats_table(outgroup="rothschildianus",set="m",chrom="X",pops=FLNC,popStr="FLNC"),
  
  stats_table(outgroup="rothschildianus",set="m",chrom="Y",pops=TXNC,popStr="TXNC"),
  stats_table(outgroup="rothschildianus",set="m",chrom="A",pops=TXNC,popStr="TXNC"),
  stats_table(outgroup="rothschildianus",set="m",chrom="H",pops=TXNC,popStr="TXNC"),
  stats_table(outgroup="rothschildianus",set="m",chrom="X",pops=TXNC,popStr="TXNC"),
  
  stats_table(outgroup="rothschildianus",set="m",chrom="Y",pops=TXFL,popStr="TXFL"),
  stats_table(outgroup="rothschildianus",set="m",chrom="A",pops=TXFL,popStr="TXFL"),
  stats_table(outgroup="rothschildianus",set="m",chrom="H",pops=TXFL,popStr="TXFL"),
  stats_table(outgroup="rothschildianus",set="m",chrom="X",pops=TXFL,popStr="TXFL")
  
 # stats_table(outgroup="rothschildianus",set="m",chrom="XY",pops=phase,popStr="phase")
  
  )
  , stringsAsFactors = FALSE)


all_data_var <- data.frame(
  rbind(
    stats_var(outgroup="bucephalophorus",set="sub4",chrom="X",pops=pop,popStr="pop",ksyn = 0.4),
    stats_var(outgroup="rothschildianus",set="sub4",chrom="X",pops=pop,popStr="pop",ksyn = 0.3),
    #Y Male
    stats_var(outgroup="bucephalophorus",set="sub4",chrom="Y",pops=pop,popStr="pop",ksyn = 0.4),
    stats_var(outgroup="rothschildianus",set="sub4",chrom="Y",pops=pop,popStr="pop",ksyn = 0.3),
    #H Male
    stats_var(outgroup="rothschildianus",set="m",chrom="H",pops=pop,popStr="pop"),
    stats_var(outgroup="bucephalophorus",set="m",chrom="H",pops=pop,popStr="pop"),
    #A Male
    stats_var(outgroup="rothschildianus",set="m",chrom="A",pops=pop,popStr="pop")
    
  )
  , stringsAsFactors = FALSE) 

chromSet <- c("A","H","X","Y")

####theta####
all_data_theta <- all_data[all_data$var == "theta" & all_data$cod == "syn"  ,]

title_thetasyn <- expression(paste(theta, ""[syn]))

ggplot(all_data_theta, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_thetasyn) +
  #theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid( pop ~ .) +
 # scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
 #   '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  ))

####pi####
title_pisyn <- expression(paste(pi, ""[syn]))
title_pi <- expression(paste(pi))

all_data_pi <- all_data[ all_data$sex == "sub4" & all_data$var == "pi" & all_data$cod == "syn" &  (all_data$pop == "XY" | all_data$pop == "XYY") ,]

#pi <- 
ggplot(all_data_pi, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_pisyn) +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(  pop ~ outgroup) +
 # facet_grid( . ~ Pop ) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',#yellow #Y
 #   '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  ))

#pop subsets (NCFL)
all_data_pi <- all_data[all_data$sex == "m" & all_data$var == "pi" & all_data$cod == "syn" & all_data$pop != "R.hastatulus" & all_data$outgroup == "rothschildianus",]
all_data_pi$Pop =factor(all_data_pi$pop, c("XY","XYY", "NC","FL"))

ggplot(all_data_pi, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_pisyn) +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(  pop ~ outgroup) +
  # facet_grid( . ~ Pop ) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',#yellow #Y
    #   '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  ))

all_data_pi <- all_data[ all_data$var == "pi" & (all_data$pop == "XY" | all_data$pop == "XYY") & all_data$outgroup == "rothschildianus",]
all_data_pi$Chrom <- factor(all_data_pi$chrom, c("A","H","X","Y"))
  
ggplot(all_data_pi, aes(x=cod, y=value)) + guides(fill = FALSE) +
    geom_bar(position=position_dodge(), stat="identity" ) +
    geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
    theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_pi) +
    theme(strip.background =element_rect(fill="white")) +
    facet_grid(  pop ~ Chrom ) 
  
####TajD####
all_data_tajD <- all_data[all_data$var == "tajD" & all_data$cod == "syn" & all_data$outgroup == "rothschildianus"& (all_data$pop == "XY" | all_data$pop == "XYY") ,]

#tajD <- 
ggplot(all_data_tajD, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) + 
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Tajima's D") +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid( . ~ pop) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  
#  scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
 #   '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  ))

multiplot(pi,tajD)

####FST####
title_fst <- expression(paste("F", ""[ST]))

all_data_fst <- all_data[ all_data$pop == "R.hastatulus" & all_data$var == "Fst" & all_data$cod == "syn"  ,]

#fst_plot <-
ggplot(all_data_fst, aes(x=chrom, y=value, fill=chrom)) +
 guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_fst) +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(. ~ outgroup) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
  #  '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  )) #+
  
  annotate(geom="text", x = "A", y = 0.1, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "H", y = 0.1, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "X", y = 0.15, label = "b", parse = TRUE, size=10) +
  annotate(geom="text", x = "Y", y = 0.45, label = "c", parse = TRUE, size=10) 

t.testloop(chromSet=chromSet,var="Fst",outgroup="rothschildianus",cod="syn",pop="R.hastatulus")

fst_data_var <- all_data_var[ all_data_var$var == "Fst", ]  
fst_anova <- aov(value ~ chrom, data=fst_data_var)
summary(fst_anova) 

FstVar <- all_data_var[all_data_var$var == "Fst" &  all_data_var$sex == "sub4" &  all_data_var$pop == "R.hastatulus" & all_data_var$cod == "syn"  ,]

#fstVar_plot <- 
ggplot(FstVar, aes(x=chrom, y=value, fill=chrom)) + 
  geom_violin()+
  geom_boxplot(width=0.1)+
  theme_bw()  + theme_bw(base_size = 30) + labs(x="",y=title_fst) +
  theme(strip.background =element_rect(fill="white")) + #ylim(0,0.01) +
  theme(strip.text = element_text(face = "italic")) +
  facet_grid(. ~ outgroup) +
  scale_x_discrete(limits=c("A","H","X","Y")) + guides(fill=FALSE) + 
  scale_fill_manual(values=c( 
    '#00A550' ,#green #A
    '#ffb14e', #orange #H
    '#00ADEF', #Blue #X
    '#FFF100'  #yellow #Y
    #   '#1B75BB', #purple-y #N
  ))

sub_fst <- all_data[all_data$sex == "m" & (all_data$pop == "FLNC" | all_data$pop == "TXNC" | all_data$pop == "TXFL" ) & all_data$outgroup == "rothschildianus" & all_data$var == "Fst" & all_data$cod == "syn"  ,]

#fst_sub_plot <-
  ggplot(sub_fst, aes(x=chrom, y=value, fill=chrom)) +
  guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_fst) +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  theme(strip.background =element_rect(fill="white")) +
   facet_grid(. ~ pop) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
    #  '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  )) 
  
####dxy####
title_dxy <- expression(paste("d", ""[XY]))
dxy <- all_data[all_data$var == "dxy" & all_data$sex == "m" & all_data$outgroup == "rothschildianus" & all_data$pop == "R.hastatulus" ,]

#dxy_plot <- 
ggplot(dxy, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  #geom_point(position=position_dodge(), stat="identity" ) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_dxy) +
  theme(strip.background =element_rect(fill="white")) +
 # facet_grid(. ~ Sex) + 
  scale_x_discrete(limits=c("A","H","X","Y")) +
  
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
  #  '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  )) #+
  annotate(geom="text", x = "A", y = 0.003, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "H", y = 0.003, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "X", y = 0.006, label = "b", parse = TRUE, size=10) +
  annotate(geom="text", x = "Y", y = 0.006, label = "b", parse = TRUE, size=10) 
  
  
t.testloop(chromSet=chromSet,var="dxy",outgroup="rothschildianus",cod="NA",pop="R.hastatulus")

dxyVar <- all_data_var[all_data_var$var == "dxy" &  all_data_var$outgroup == "rothschildianus" & all_data_var$pop == "R.hastatulus" ,]

#dxyVar_plot <- 
ggplot(dxyVar, aes(x=chrom, y=value, fill=chrom)) + 
    geom_violin()+
  geom_boxplot(width=0.1)+
  theme_bw()  + theme_bw(base_size = 30) + labs(x="",y=title_dxy) +
  theme(strip.background =element_rect(fill="white")) + #ylim(0,0.01) +
  theme(strip.text = element_text(face = "italic")) +
  scale_x_discrete(limits=c("A","H","X","Y")) + guides(fill=FALSE) + 
  scale_fill_manual(values=c( 
    '#00A550' ,#green #A
    '#ffb14e', #orange #H
    '#00ADEF', #Blue #X
    '#FFF100'  #yellow #Y
    #   '#1B75BB', #purple-y #N
  ))
  
multiplot(fst_plot,dxy_plot,cols=2)
multiplot(fstVar_plot,dxyVar_plot,cols=2)

sub_dxy <- all_data[all_data$sex == "m" & (all_data$pop == "FLNC" | all_data$pop == "TXNC" | all_data$pop == "TXFL" ) & all_data$outgroup == "rothschildianus" & all_data$var == "dxy"   ,]

#dxy_sub_plot <-
ggplot(sub_dxy, aes(x=chrom, y=value, fill=chrom)) +
  guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_dxy) +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(. ~ pop) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
    #  '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  )) 

multiplot(fst_sub_plot,dxy_sub_plot)

####da####
da <- all_data[all_data$var == "da" & (all_data$pop == "FLNC" | all_data$pop == "TXNC" | all_data$pop == "TXFL" ) & all_data$outgroup == "rothschildianus",]

#da <- all_data[all_data$var == "da" & all_data$pop == "R.hastatulus" & all_data$outgroup == "rothschildianus",]

title_da <- expression(paste("d", ""[A]))

#da_plot <- 
  ggplot(da, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  #geom_point(position=position_dodge(), stat="identity" ) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_da) +
  theme(strip.background =element_rect(fill="white")) +
   facet_grid(. ~ pop) + 
  #  scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
    #  '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  ))  #+
  
  annotate(geom="text", x = "A", y = 0.0005, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "H", y = 0.00075, label = "b", parse = TRUE, size=10) +
  annotate(geom="text", x = "X", y = 0.00125, label = "c", parse = TRUE, size=10) +
  annotate(geom="text", x = "Y", y = 0.0026, label = "d", parse = TRUE, size=10) 
  
  
  multiplot(fst_plot,dxy_plot,da_plot,cols=3)
  multiplot(fst_sub_plot,dxy_sub_plot,da_plot)
  
####Ks####
ks <- all_data[all_data$var == "k" & (all_data$pop == "XY" | all_data$pop == "XYY")  & all_data$cod == "syn",]

ks$Chromosome = factor(ks$chrom, levels=c("A","H","X","Y"))
#ks$Outgroup = factor(ks$outgroup, levels=c("R.rothschildianus", "R.bucephalophorus"))

#ks_plot <- 
ggplot(ks, aes(x=Chromosome, y=value, fill=Chromosome )) + 
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + guides(fill=FALSE) +
  theme_bw()  + theme_bw(base_size = 30) + labs(x="",y="Synonymous Divergence (Ks)") +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(. ~ outgroup) +
  scale_x_discrete(limits=c("A","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00A550', #green #A
#    '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00ADEF', #Blue #X
    '#FFF100'  #yellow #Y
  ))


####dnds####
evo_rate <-  all_data[all_data$var == "kaks" & all_data$pop == "R.hastatulus",]
  
evo_rate$outgroup[evo_rate$outgroup == "rothschildianus"] <- "R.rothschildianus"
evo_rate$outgroup[evo_rate$outgroup == "bucephalophorus"] <- "R.bucephalophorus"

#rate_plot <-
ggplot(evo_rate, aes(x=chrom, y=value, color=chrom)) + guides(color=FALSE) +
  geom_point(position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Divergence (dn/ds)") +
  theme(strip.background =element_rect(fill="white")) +
  scale_x_discrete(limits=c("A","H","X","Y")) + 
  facet_grid(. ~ outgroup) +
  scale_color_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
 #   '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  ))

kaksvar <- all_data_var[all_data_var$var == "kaks" & all_data_var$pop == "R.hastatulus" ,]
kaksvar$outgroup[kaksvar$outgroup == "rothschildianus"] <- "R.rothschildianus"
kaksvar$outgroup[kaksvar$outgroup == "bucephalophorus"] <- "R.bucephalophorus"

#rate_plot <-
ggplot(kaksvar, aes(x=chrom, y=value, fill=chrom)) + 
  geom_abline(intercept = 1,slope=0,alpha=0.2) +
  geom_violin()+
  geom_boxplot(width=0.1)+
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Divergence (dn/ds)") +
  theme(strip.background =element_rect(fill="white")) + ylim(0,10) +
  theme(strip.text = element_text(face = "italic")) +
  scale_x_discrete(limits=c("A","H","X","Y")) + guides(fill=FALSE) + #coord_trans( y = "log2") +
  facet_grid(. ~ outgroup) +
  scale_fill_manual(values=c( 
    '#00A550' ,#green #A
    '#ffb14e', #orange #H
    '#00ADEF', #Blue #X
    '#FFF100'  #yellow #Y
    #   '#1B75BB', #purple-y #N
    
    
  )) #+ 
  
  annotate(geom="text", x = "A", y = 2, label = '"a"', parse = TRUE, size=10) +
  annotate(geom="text", x = "H", y = 2, label = '"b"', parse = TRUE, size=10) +
  annotate(geom="text", x = "X", y = 2, label = '"a"', parse = TRUE, size=10) +
  annotate(geom="text", x = "Y", y = 2, label = '"a"', parse = TRUE, size=10) 

t.testloop(chromSet=chromSet,var="kaks",pop="R.hastatulus",cod="NA",outgroup="bucephalophorus")
t.testloop(chromSet=chromSet,var="kaks",pop="R.hastatulus",cod="NA",outgroup="rothschildianus")

all_data_k <- all_data[ all_data$var == "k" & (all_data$pop == "XY" | all_data$pop == "XYY") & all_data$outgroup == "rothschildianus",]
all_data_k$Chrom <- factor(all_data_pi$chrom, c("A","H","X","Y"))

ggplot(all_data_k, aes(x=cod, y=value)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="k") +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(  pop ~ Chrom ) 

t.testloop(chromSet=chromSet,var="k",pop="XY",cod="rep",outgroup="rothschildianus")

  
####mk####

mk <- all_data[all_data$var == "mk" & all_data$pop == "R.hastatulus" & all_data$sex == "m"  ,]
mk$outgroup[mk$outgroup == "rothschildianus"] <- "R.rothschildianus"
mk$outgroup[mk$outgroup == "bucephalophorus"] <- "R.bucephalophorus"

#mk$Chromosome = factor(mk$chrom, levels=c('A','XY','N','H'))
#mk$Chromosome = factor(mk$chrom, levels=c("A","N","H","X","Y"))

mk$Outgroup =factor(mk$outgroup, c("R.rothschildianus", "R.bucephalophorus"))

#adapt_plot <-
ggplot(mk, aes(x=chrom, y=value, color=chrom)) + 
  geom_abline(intercept = 1,slope=0,alpha=0.2) +
  geom_point(position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Adaptation (MK)") +
  theme(strip.background =element_rect(fill="white")) + #ylim(0.5,2.1) +
 # scale_x_discrete(limits=c("A","N","H","X","Y")) + guides(color=FALSE) +
  scale_x_discrete(limits=c("A","H","X","Y")) + guides(color=FALSE) + coord_trans( y = "log2")+
   theme(strip.text = element_text(face = "italic")) +
  facet_grid(. ~ Outgroup) +
  scale_color_manual(values=c( 
    '#00ADEF', #Blue #X
    '#FFF100',  #yellow #Y
 #   '#1B75BB', #purple-y #N
    '#ffb14e', #orange #H
    '#00A550' #green #A
  )) + 
  
  annotate(geom="text", x = "H", y = 2.1, label = '"*"', parse = TRUE, size=10) 

mkvar <- all_data_var[all_data_var$var == "mk" & all_data_var$pop != "R.hastatulus"  ,]
mkvar$outgroup[mkvar$outgroup == "rothschildianus"] <- "R.rothschildianus"
mkvar$outgroup[mkvar$outgroup == "bucephalophorus"] <- "R.bucephalophorus"

#adapt_plot <-
ggplot(mkvar[mkvar$locus %nin% bad_dnds_div,], aes(x=chrom, y=value, fill=chrom)) + 
#  ggplot(mkvar, aes(x=chrom, y=value, fill=chrom)) + 
  
    geom_abline(intercept = 1,slope=0,alpha=0.2) +
  geom_violin()+
  geom_boxplot(width=0.1)+
    theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Adaptation (MK)") +
  theme(strip.background =element_rect(fill="white")) + ylim(0.001,100) +
   theme(strip.text = element_text(face = "italic")) +
  scale_x_discrete(limits=c("A","H","X","Y")) + guides(fill=FALSE) + coord_trans( y = "log10")+
  facet_grid( pop ~ outgroup) +
  scale_fill_manual(values=c( 
    '#00A550' ,#green #A
    '#ffb14e', #orange #H
     '#00ADEF', #Blue #X
    '#FFF100'  #yellow #Y
    #   '#1B75BB', #purple-y #N
    
    
  )) #+ 
  
  annotate(geom="text", x = "A", y = 4, label = '"a"', parse = TRUE, size=10) +
  annotate(geom="text", x = "H", y = 4, label = '"ab"', parse = TRUE, size=10) +
  annotate(geom="text", x = "X", y = 4, label = '"b"', parse = TRUE, size=10) +
  annotate(geom="text", x = "Y", y = 4, label = '"c"', parse = TRUE, size=10) 
  
  t.testloop(chromSet=chromSet,var="mk",pop="R.hastatulus",cod="NA",outgroup="bucephalophorus",eLoc=bad_dnds_div)
  t.testloop(chromSet=chromSet,var="mk",pop="R.hastatulus",cod="NA",outgroup="rothschildianus",eLoc=bad_dnds_div)
  
  t.testloop(chromSet=chromSet,var="mk",pop="XY",cod="NA",outgroup="bucephalophorus",eLoc=bad_dnds_div)
  t.testloop(chromSet=chromSet,var="mk",pop="XY",cod="NA",outgroup="rothschildianus",eLoc=bad_dnds_div)
  
  t.testloop(chromSet=chromSet,var="mk",pop="XYY",cod="NA",outgroup="bucephalophorus",eLoc=bad_dnds_div)
  t.testloop(chromSet=chromSet,var="mk",pop="XYY",cod="NA",outgroup="rothschildianus",eLoc=bad_dnds_div)
  


multiplot(rate_plot,adapt_plot)
####Q####
##Qpi##

Q_data_pi <- data.frame( cbind(
    c("XY","XYY"),
    c(
    all_data_pi$value[all_data_pi$chrom == "X" & all_data_pi$pop == "XY" & all_data_pi$outgroup == "rothschildianus"  ] / all_data_pi$value[all_data_pi$chrom == "A" & all_data_pi$pop == "XY" & all_data_pi$outgroup == "rothschildianus" ],
    all_data_pi$value[all_data_pi$chrom == "X" & all_data_pi$pop == "XYY" & all_data_pi$outgroup == "rothschildianus"   ] / all_data_pi$value[all_data_pi$chrom == "A" & all_data_pi$pop == "XYY" & all_data_pi$outgroup == "rothschildianus"  ]
    ),
    c("pi","pi")
  ) 
)
names(Q_data_pi) <- c("Pop","Q","type")



options(digits = 4)
Q_data_pi$Q<- type.convert(Q_data_pi$Q, na.strings = "NA", as.is = FALSE, dec = ".",
             numerals =  "no.loss")

##Qfst##

Q_data_fst <-
  data.frame(cbind(
      "R.hastatulus",
      log((1- (2*all_data_fst$value[all_data_fst$chrom == "X" & all_data_fst$outgroup == "rothschildianus" ]))) / log((1-(2*all_data_fst$value[all_data_fst$chrom == "A" & all_data_fst$outgroup == "rothschildianus" ]))),
      "Fst"
      )
    )

names(Q_data_fst) <- c("Pop","Q","type")
Q_data_fst$Q<- type.convert(Q_data_fst$Q, na.strings = "NA", as.is = FALSE, dec = ".", numerals =  "no.loss")

##Qdxy##

Q_data_dxy <-
  data.frame(cbind(
    "R.hastatulus",
    dxy$value[dxy$chrom == "X" & dxy$outgroup == "rothschildianus" ] /  dxy$value[dxy$chrom == "A" & dxy$outgroup == "rothschildianus" ],
    "Dxy"
  )
  )

names(Q_data_dxy) <- c("Pop","Q","type")
Q_data_dxy$Q<- type.convert(Q_data_dxy$Q, na.strings = "NA", as.is = FALSE, dec = ".", numerals =  "no.loss")


Q<-rbind(Q_data_fst,Q_data_pi,Q_data_dxy)

#title_qfst <- expression(paste("Q", ""[Fst]))

ggplot(Q, aes(x=type, y=Q,fill=Pop)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Q") 



####MS####
#import S for all
#within <- fread('m_rothschildianus_summarystats_pop_Y.txt')
#hast_S <- data.frame( cbind( within$locus, within$pop0_sites_syn ))
#hast_S <- hast_S[hast_S$X2 != 0,]
#write.csv(hast_S, file = "hast_synsites.txt" )
#write(hast_S, file = "hast_synsites.txt", append = FALSE, sep = "\n")




##pi##
ms_pi <- rbind(
               ms_stat(chrom="A",var="pi_tot",sitemean=265.798,filename = "ancMig.parsed.ms"),
               cbind(all_data_pi[all_data_pi$pop == "R.hastatulus",-c(2,3)],state="obs")
)

#pi_plot_oe <-
ggplot(ms_pi, aes(x=chrom, y=value, fill=state)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_pisyn) #+
  #scale_x_discrete(limits=c("A","H","X","Y")) 
  scale_x_discrete(limits=c("A","X")) 

##fst##
ms_fst <- rbind(#ms_stat(chrom="X",var="fst"),
                ms_stat(chrom="A",var="fst",filename = "ancMig.parsed.ms"),
                cbind(all_data_fst[all_data_fst$pop == "R.hastatulus",-c(2,3)],state="obs")
)

#fst_plot_oe <- 
  ggplot(ms_fst, aes(x=chrom, y=value, fill=state)) +
  #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Fst") #+
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  #scale_x_discrete(limits=c("A","H","X","Y")) 
  scale_x_discrete(limits=c("A","X")) 


##dxy##

ms_dxy <- rbind(
  #ms_stat(chrom="X",var="dxy",sitemean=2.50881),
  ms_stat(chrom="A",var="dxy",sitemean=2.50881,filename = "ancMig.parsed.ms"),
  cbind(kxydxy[kxydxy$pop == "R.hastatulus",-c(2,3)],"state"="obs")
)

#dxy_plot_oe <- 
  ggplot(ms_dxy, aes(x=chrom, y=value, fill=state)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Dxy") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  #scale_x_discrete(limits=c("A","H","X","Y")) 
  scale_x_discrete(limits=c("A","X")) 


####Extra Calculations####
##synsite##
all_data[all_data$var == "sites" & all_data$cod == "syn" 
         & all_data$pop == "R.hastatulus" & all_data$chrom == "A", ]

#NE 

0.00629 / (4 * 7.5e-9) # Ne ancestral

(0.00629 / (4 * 7.5e-9) ) * 3.4 
(0.00629 / (4 * 7.5e-9) ) * 0.079


####X2Y#####

X2Ydata<-cbind(fread('m_rothschildianus_summarystats_phase_XY.txt'),fread('m_rothschildianus_interpop_phase_XY.txt'))
X2Ydata <- X2Ydata[,-c(51)]

X <- ggplot(X2Ydata,aes(x=pop0_d_syn,y=pop1_kaks_NA)) + geom_point() + ylim(0,2)

Y <- ggplot(X2Ydata,aes(x=pop0_d_syn,y=pop2_kaks_NA)) + geom_point() + ylim(0,2)

multiplot(X,Y)

##cytotype specific mk####


rxyyh <- cbind(fread('xyy_rothschildianus_summarystats_pop_h.txt'),"pop"="xyy","outgroup"="rothschildianus","chrom"="h")
rxyya <- cbind(fread('xyy_rothschildianus_summarystats_pop_a.txt'),"pop"="xyy","outgroup"="rothschildianus","chrom"="a")

rxyh <- cbind(fread('xy_rothschildianus_summarystats_pop_h.txt'),"pop"="xy","outgroup"="rothschildianus","chrom"="h")
rxya <- cbind(fread('xy_rothschildianus_summarystats_pop_a.txt'),"pop"="xy","outgroup"="rothschildianus","chrom"="a")

bxyyh <- cbind(fread('xyy_bucephalophorus_summarystats_pop_h.txt'),"pop"="xyy","outgroup"="bucephalophorus","chrom"="h")
bxyya <- cbind(fread('xyy_bucephalophorus_summarystats_pop_a.txt'),"pop"="xyy","outgroup"="bucephalophorus","chrom"="a")

bxyh <- cbind(fread('xy_bucephalophorus_summarystats_pop_h.txt'),"pop"="xy","outgroup"="bucephalophorus","chrom"="h")
bxya <- cbind(fread('xy_bucephalophorus_summarystats_pop_a.txt'),"pop"="xy","outgroup"="bucephalophorus","chrom"="a")

rxyyh <- rxyyh[rxyyh$pop0_k_syn < 0.3,]
rxyya <- rxyya[rxyya$pop0_k_syn < 0.3,]
rxyh <- rxyh[rxyh$pop0_k_syn < 0.3,]
rxya <- rxya[rxya$pop0_k_syn < 0.3,]

bxyyh <- bxyyh[bxyyh$pop0_k_syn < 0.4,]
bxyya <- bxyya[bxyya$pop0_k_syn < 0.4,]
bxyh <- bxyh[bxyh$pop0_k_syn < 0.4,]
bxya <- bxya[bxya$pop0_k_syn < 0.4,]

filtered <- rbind(
rxyyh[rxyyh$pop0_k_syn < 0.3,],
rxyya[rxyya$pop0_k_syn < 0.3,],
rxyh[rxyh$pop0_k_syn < 0.3,],
rxya[rxya$pop0_k_syn < 0.3,],
bxyyh[bxyyh$pop0_k_syn < 0.4,],
bxyya[bxyya$pop0_k_syn < 0.4,],
bxyh[bxyh$pop0_k_syn < 0.4,],
bxya[bxya$pop0_k_syn < 0.4,],fill=TRUE
)

ggplot(filtered, aes(x=chrom, y=pop0_kaks_NA, fill=chrom)) + 
  #  ggplot(mkvar, aes(x=chrom, y=value, fill=chrom)) + 
  geom_abline(intercept = 1,slope=0,alpha=0.2) +
  geom_violin()+
  geom_boxplot(width=0.1)+
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Substitution (dnds)") +
  theme(strip.background =element_rect(fill="white")) + 
  theme(strip.text = element_text(face = "italic")) +
  #scale_x_discrete(limits=c("A","H","X","Y")) + 
  guides(fill=FALSE) + coord_trans( y = "log2")+
  facet_grid( pop ~ outgroup) 


ggplot(filtered, aes(x=chrom, y=pop0_mk_NA, fill=chrom)) + 
  #  ggplot(mkvar, aes(x=chrom, y=value, fill=chrom)) + 
  geom_abline(intercept = 1,slope=0,alpha=0.2) +
  geom_violin()+
  geom_boxplot(width=0.1)+
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Adaptation (MK)") +
  theme(strip.background =element_rect(fill="white")) + 
  theme(strip.text = element_text(face = "italic")) +
  #scale_x_discrete(limits=c("A","H","X","Y")) + 
  guides(fill=FALSE) + coord_trans( y = "log2")+
  facet_grid( pop ~ outgroup) 

t.test(
  log(filtered$pop0_kaks_NA[filtered$pop == "xyy" & filtered$outgroup == "bucephalophorus" & filtered$chrom == "h" ]),
  log(filtered$pop0_kaks_NA[filtered$pop == "xyy" & filtered$outgroup == "bucephalophorus" & filtered$chrom == "a" ])
       )

t.test(
  log(filtered$pop0_kaks_NA[filtered$pop == "xy" & filtered$outgroup == "bucephalophorus" & filtered$chrom == "h" ]),
  log(filtered$pop0_kaks_NA[filtered$pop == "xy" & filtered$outgroup == "bucephalophorus" & filtered$chrom == "a" ])
)

t.test(
  log(filtered$pop0_kaks_NA[filtered$pop == "xyy" & filtered$outgroup == "rothschildianus" & filtered$chrom == "h" ]),
  log(filtered$pop0_kaks_NA[filtered$pop == "xyy" & filtered$outgroup == "rothschildianus" & filtered$chrom == "a" ])
)

t.test(
  log(filtered$pop0_kaks_NA[filtered$pop == "xy" & filtered$outgroup == "rothschildianus" & filtered$chrom == "h" ]),
  log(filtered$pop0_kaks_NA[filtered$pop == "xy" & filtered$outgroup == "rothschildianus" & filtered$chrom == "a" ])
)
