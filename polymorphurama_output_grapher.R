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

stats_table <- function(outgroup=NULL,set=NULL,chrom=NULL,pops=NULL,subsetList=NULL,popStr='pop'){ 
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

    filename_wIn <- paste(set,"_",outgroup,"_summarystats_",popStr,"_",chrom,".txt",sep="")
    filename_btw <- paste(set,"_",outgroup,"_interpop_",popStr,"_",chrom,".txt",sep="")
  in_read <- fread(filename_wIn)
  #in_read$pop0_pi_syn_w <- in_read$pop0_pi_syn * in_read$pop0_sites_syn
  seqMax <- max(in_read$pop0_seqs_NA[!is.na(in_read$pop0_seqs_NA)])
  coverage_list <- in_read$locus[in_read$pop0_seqs_NA >= seqMax  & in_read$pop0_sites_syn >= 50]
  in_read_cov <- in_read[in_read$locus %in% coverage_list ]
  in_bet <- fread(filename_btw)
  in_bet_cov <- in_bet[in_bet$locus %in% coverage_list ]
  if (length(subsetList) >0 ){
    in_read_cov <- in_read[in_read$locus %in% subsetList ]
    in_bet_cov <- in_bet_cov[in_bet_cov$locus %in% subsetList ]
  }
  stats_wIn <- summaryStats(in_mat=in_read_cov ,pops=pops,chrom=chrom)
  stats_btw <- summaryStatsInter(inter=in_bet_cov,pops=pops,chrom=chrom)
  stats_bind <- cbind(rbind(stats_wIn,stats_btw),outgroup=outgroup,stringsAsFactors=FALSE)
  return(stats_bind)
}
stats_var <- function(outgroup=NULL,set=NULL,chrom=NULL,pops=NULL,popStr='pop'){ 

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
  filename_wIn <- paste(set,"_",outgroup,"_summarystats_",popStr,"_",chrom,".txt",sep="")
  filename_btw <- paste(set,"_",outgroup,"_interpop_",popStr,"_",chrom,".txt",sep="")
  
  in_read <- fread(filename_wIn)
  #in_read$pop0_pi_syn_w <- in_read$pop0_pi_syn * in_read$pop0_sites_syn
  seqMax <- max(in_read$pop0_seqs_NA[!is.na(in_read$pop0_seqs_NA)])
  coverage_list <- in_read$locus[in_read$pop0_seqs_NA == seqMax & in_read$pop0_sites_syn >= 100]
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


####import####

pop <- c("R.hastatulus","XY","XYY")
FLNC <- c("XYY","FL","NC")

all_data <- data.frame(
  rbind(
    stats_table(outgroup="bucephalophorus",set="rna",chrom="H",pops=pop,popStr="pop"),
    stats_table(outgroup="rothschildianus",set="rna",chrom="H",pops=pop,popStr="pop"),
    stats_table(outgroup="bucephalophorus",set="rna",chrom="XY",pops=pop,popStr="pop"),
    stats_table(outgroup="rothschildianus",set="rna",chrom="XY",pops=pop,popStr="pop"),
    stats_table(outgroup="bucephalophorus",set="rna",chrom="N",pops=pop,popStr="pop"),
    stats_table(outgroup="rothschildianus",set="rna",chrom="N",pops=pop,popStr="pop"),
    stats_table(outgroup="rothschildianus",set="rna",chrom="A",pops=pop,popStr="pop")
  )
, stringsAsFactors = FALSE)

all_data_var <- data.frame(
  rbind(
    stats_var(outgroup="rothschildianus",set="josh",chrom="A",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="josh",chrom="XY",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="josh",chrom="N",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="josh",chrom="H",pops=pop,popStr="pop")
  )
, stringsAsFactors = FALSE)

####theata####
all_data_theta <- all_data[all_data$var == "theta" & all_data$cod == "syn"  ,]

title_thetasyn <- expression(paste(theta, ""[syn]))

ggplot(all_data_theta, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_thetasyn) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  facet_grid(. ~ pop, scales = "free") +
  scale_x_discrete(limits=c("A","XY","H","N")) #+
scale_fill_manual(values=c( 
  '#00ADEF', #Blue
  '#FFF100',  #yellow
  '#00A550', #green
  '#1B75BB', #purple-y
  '#8B4BD8'
))

####pi####

all_data_pi <- all_data[all_data$var == "pi" & all_data$cod == "syn"  ,]

title_pisyn <- expression(paste(pi, ""[syn]))

ggplot(all_data_pi, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_pisyn) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  facet_grid(. ~ pop, scales = "free") +
  scale_x_discrete(limits=c("A","XY","H","N")) #+
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#1B75BB', #purple-y
    '#8B4BD8'
  ))


####TajD####

all_data_tajD <- all_data[all_data$var == "tajD" & all_data$cod == "syn"  ,]

ggplot(all_data_tajD, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) + 
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Tajima's D") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  facet_grid(. ~ pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  scale_x_discrete(limits=c("A","H","N","XY")) +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text = element_text(face = "italic")) #+
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#1B75BB', #purple-y
    '#8B4BD8'
  ))


####FST####

all_data_fst <- all_data[all_data$var == "Fst" & all_data$cod == "syn"  ,]

#fst_plot <-
ggplot(all_data_fst, aes(x=chrom, y=value, fill=chrom)) +
 guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Fst") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  scale_x_discrete(limits=c("A","XY","H","N")) #+
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550' #green
   #'#1B75BB', #purple-y
   # '#8B4BD8'
  )) +
  
  annotate(geom="text", x = "A", y = 0.12, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "X", y = 0.155, label = "b", parse = TRUE, size=10) +
  annotate(geom="text", x = "Y", y = 0.51, label = "c", parse = TRUE, size=10) 


####dxy####

dxy <- all_data[all_data$var == "dxy" ,]
#dxy <- dxy[c(1,2,4),]

#dxy_plot <- 
ggplot(dxy, aes(x=chrom, y=value, fill=chrom)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Dxy") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  scale_x_discrete(limits=c("A","XY","N","H")) #+
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550' #green

  )) +
  annotate(geom="text", x = "A", y = 0.014, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "X", y = 0.006, label = "b", parse = TRUE, size=10) +
  annotate(geom="text", x = "Y", y = 0.007, label = "b", parse = TRUE, size=10) 
  
multiplot(fst_plot,dxy_plot,cols=2)


####kxy####
kxydxy <- rbind(
all_data[all_data$var == "kxy" & all_data$pop == "R.hastatulus",],
all_data[ all_data$var == "dxy",]
)
kxydxy$outgroup[kxydxy$var == "dxy"] <- "R.hastatulus"
kxydxy$outgroup[kxydxy$outgroup == "rothschildianus"] <- "R.rothschildianus"

#levels(kxydxy$outgroup) = c("R.rothschildianus", "R.hastatulus","R.bucephalophorus")
#levels(kxydxy$outgroup) = c("R.rothschildianus", "R.hastatulus")


#div_plot <- 
ggplot(kxydxy, aes(x=outgroup, y=value, color=chrom
                      # , group=1
                       )) + #guides(fill = FALSE) +
  geom_point(position=position_dodge(.9)) +
  
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Average Divergence") +
  scale_x_discrete(limits=c("R.hastatulus","R.rothschildianus","R.bucephalophorus")) +
   theme(axis.text.x = element_text(face = "italic"))


####dnds####
evo_rate <- rbind(
  all_data[all_data$var == "kaks" & all_data$pop == "R.hastatulus",],
  all_data[all_data$var == "dnds" ,]
  
)

#levels(all_data_dxy$outgroup) = c("rothschildianus", "hastatulus","bucephalophorus")
evo_rate$outgroup[evo_rate$var == "dnds"] <- "R.hastatulus"
evo_rate$outgroup[evo_rate$outgroup == "rothschildianus"] <- "R.rothschildianus"

#rate_plot <-
ggplot(evo_rate, aes(x=outgroup, y=value, color=chrom
                         # , group=1
)) + #guides(fill = FALSE) +
  geom_point(position=position_dodge(.9)) +
  
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Divergence") +
  scale_x_discrete(limits=c("R.hastatulus","R.rothschildianus","R.bucephalophorus")) +
  theme(axis.text.x = element_text(face = "italic"))


####mk####

mk <- all_data[all_data$var == "mk" & all_data$pop == "R.hastatulus",]
mk$outgroup[mk$outgroup == "rothschildianus"] <- "R.rothschildianus"

mk$Chromosome = factor(mk$chrom, levels=c('A','XY','N','H'))

#adapt_plot <-
ggplot(mk, aes(x=outgroup, y=value, color=Chromosome
                     # , group=1
)) + #guides(fill = FALSE) +
  geom_point(position=position_dodge(.9)) +

  geom_abline(intercept = 1,slope=0,alpha=0.5) +
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Adaptation") +
  scale_x_discrete(limits=c("R.rothschildianus","R.bucephalophorus")) +
  theme(axis.text.x = element_text(face = "italic")) + coord_trans(y = "log2")

####MS####



#filename <- paste('ms_',chrom,'_stat.txt',sep="")


##pi##
ms_pi <- rbind(
               ms_stat(chrom="A",var="pi_tot",sitemean=265.798,filename = "secMig.parsed.ms"),
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
                ms_stat(chrom="A",var="fst",filename = "secMig.parsed.ms"),
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

t.test(
  all_data_var$value[all_data_var$chrom == "A" & all_data_var$var == "Fst" ] ,  
  all_data_var$value[all_data_var$chrom == "X" & all_data_var$var == "Fst" ])

t.test(
  all_data_var$value[all_data_var$chrom == "N" & all_data_var$var == "Fst" ] ,  
  all_data_var$value[all_data_var$chrom == "Y" & all_data_var$var == "Fst" ])

t.test(
  all_data_var$value[all_data_var$chrom == "X" & all_data_var$var == "Fst" ] ,  
  all_data_var$value[all_data_var$chrom == "Y" & all_data_var$var == "Fst" ])

#Nem = ¼(fst - 1)  
(1/0.1165208 - 1)/4

fst_data_var <- all_data_var[ all_data_var$var == "Fst", ]  
fst_anova <- aov(value ~ chrom, data=fst_data_var)
summary(fst_anova) 

##dxy##


ms_dxy <- rbind(
  ms_stat(chrom="X",var="dxy",sitemean=2.50881),
  ms_stat(chrom="A",var="dxy",sitemean=2.50881),
  cbind(dxy[,-c(2,3)],"state"="obs")
)

dxy_plot_oe <- 
  ggplot(ms_dxy, aes(x=chrom, y=value, fill=state)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Dxy") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  #scale_x_discrete(limits=c("A","H","X","Y")) 
  scale_x_discrete(limits=c("A","X")) 

multiplot(pi_plot_oe,fst_plot_oe,dxy_plot_oe,cols=3)

t.test(
  all_data_var$value[all_data_var$chrom == "A" & all_data_var$var == "dxy" ] ,  
  all_data_var$value[all_data_var$chrom == "X" & all_data_var$var == "dxy" ])

t.test(
  all_data_var$value[all_data_var$chrom == "Y" & all_data_var$var == "dxy" ] ,  
  all_data_var$value[all_data_var$chrom == "X" & all_data_var$var == "dxy" ])


dxy_data_var <- all_data_var[ all_data_var$var == "dxy", ]  
dxy_anova <- aov(value ~ chrom, data=dxy_data_var)
summary(dxy_anova) 

#Dxy = 2μt
0.014862532/(2*7.5e-9)
(115000+35000)-0.014862532/(2*7.5e-9)

####Extra Calculations####
##fmaleX!=maleX##
maleX <- data.frame(stats_table(outgroup="rothschildianus",set="XYphased",chrom="X",pops=pop))
femaleX <- data.frame(stats_table(outgroup="rothschildianus",set="joshFem",chrom="X",pops=pop))

maleX[maleX$var == "pi",]
femaleX[femaleX$var == "pi",]

##synsite##
all_data[all_data$var == "sites" & all_data$cod == "syn" 
         & all_data$pop == "R.hastatulus" & all_data$chrom == "A", ]

#Dxy = 2μt
0.269506300/(2*7.5e-9)



