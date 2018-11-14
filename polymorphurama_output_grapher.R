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
  
  ##Calculate summaries function
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
  coverage_list <- in_read$locus[in_read$pop0_seqs_NA >= (seqMax * 0.8)  & in_read$pop0_sites_syn >= 50]
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
mpop <- c("R.hastatulusMale","XYMale","XYYMale")
fpop <- c("R.hastatulusFemale","XYFemale","XYYFemale")

all_data <- data.frame(
  rbind(
    #Xchrom Male 
    stats_table(outgroup="bucephalophorus",set="m",chrom="X",pops=pop,popStr="pop"),
    stats_table(outgroup="rothschildianus",set="m",chrom="X",pops=pop,popStr="pop"),
    #Ychrom Male
    stats_table(outgroup="bucephalophorus",set="m",chrom="Y",pops=pop,popStr="pop"),
    stats_table(outgroup="rothschildianus",set="m",chrom="Y",pops=pop,popStr="pop"),
    #A Male
    stats_table(outgroup="rothschildianus",set="m",chrom="A",pops=pop,popStr="pop"),
    #H Male
    stats_table(outgroup="rothschildianus",set="m",chrom="H",pops=pop,popStr="pop"),
    stats_table(outgroup="bucephalophorus",set="m",chrom="H",pops=pop,popStr="pop"),
    #Xchrom Female
    stats_table(outgroup="bucephalophorus",set="f",chrom="X",pops=pop,popStr="pop"),
    stats_table(outgroup="rothschildianus",set="f",chrom="X",pops=pop,popStr="pop"),
    #N Female
    stats_table(outgroup="bucephalophorus",set="f",chrom="N",pops=pop,popStr="pop"),
    stats_table(outgroup="rothschildianus",set="f",chrom="N",pops=pop,popStr="pop"),
    #A Female
    stats_table(outgroup="rothschildianus",set="f",chrom="A",pops=pop,popStr="pop"),
    stats_table(outgroup="bucephalophorus",set="f",chrom="A",pops=pop,popStr="pop"),
    #H Female
    stats_table(outgroup="rothschildianus",set="f",chrom="H",pops=pop,popStr="pop"),
    stats_table(outgroup="bucephalophorus",set="f",chrom="H",pops=pop,popStr="pop")
  )
  , stringsAsFactors = FALSE)

all_data$Sex <- levels(all_data$sex, c("Male", "Female"))
#all_data$sex[all_data$sex == "f"] <- "Female"
#all_data$sex[all_data$sex == "m"] <- "Male"

all_data_var <- data.frame(
  rbind(
    stats_var(outgroup="bucephalophorus",set="rna",chrom="Xm",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="rna",chrom="Xm",pops=pop,popStr="pop"),
    stats_var(outgroup="bucephalophorus",set="rna",chrom="Ychrom",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="rna",chrom="Ychrom",pops=pop,popStr="pop"),
    stats_var(outgroup="bucephalophorus",set="XY",chrom="Xf",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="XY",chrom="Xf",pops=pop,popStr="pop"),
    stats_var(outgroup="bucephalophorus",set="N",chrom="N",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="N",chrom="N",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="A",chrom="A",pops=pop,popStr="pop"),
    stats_var(outgroup="rothschildianus",set="rna",chrom="H",pops=pop,popStr="pop"),
    stats_var(outgroup="bucephalophorus",set="rna",chrom="H",pops=pop,popStr="pop")
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
  #theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(sex ~ pop) +
  scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#ffb14e', #orange
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
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(pop ~ sex) +
  scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#ffb14e', #orange
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
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(sex ~ pop) +
  scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#ffb14e', #orange
    '#1B75BB', #purple-y
    '#8B4BD8'
  ))


####FST####

all_data_fst <- all_data[all_data$var == "Fst" & all_data$cod == "syn"  ,]
#all_data_fst <- all_data[all_data$var == "Fst" & all_data$cod == "syn" & 
#                           (all_data$chrom == "Xm" | all_data$chrom == "A" | all_data$chrom == "Ychrom") ,]


#fst_plot <-
ggplot(all_data_fst, aes(x=chrom, y=value, fill=chrom)) +
 guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Fst") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(. ~ sex) +
  scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#ffb14e', #orange
    '#1B75BB', #purple-y
    '#8B4BD8'
  ))

  
  annotate(geom="text", x = "A", y = 0.1, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "Xm", y = 0.2, label = "b", parse = TRUE, size=10) +
  annotate(geom="text", x = "Ychrom", y = 0.5, label = "c", parse = TRUE, size=10) 

  t.test(
    all_data_var$value[all_data_var$chrom == "A" & all_data_var$var == "Fst" ] ,  
    all_data_var$value[all_data_var$chrom == "Xm" & all_data_var$var == "Fst" ])
  
  t.test(
    all_data_var$value[all_data_var$chrom == "A" & all_data_var$var == "Fst" ] ,  
    all_data_var$value[all_data_var$chrom == "Ychrom" & all_data_var$var == "Fst" ])
  
  t.test(
    all_data_var$value[all_data_var$chrom == "Ychrom" & all_data_var$var == "Fst" ] ,  
    all_data_var$value[all_data_var$chrom == "Xm" & all_data_var$var == "Fst" ])
  
  fst_data_var <- all_data_var[ all_data_var$var == "Fst", ]  
  fst_anova <- aov(value ~ chrom, data=fst_data_var)
  summary(fst_anova) 
  
####dxy####

dxy <- all_data[all_data$var == "dxy" ,]
#  dxy <- all_data[all_data$var == "dxy" & 
#                    (all_data$chrom == "Xm" | all_data$chrom == "A" | all_data$chrom == "Ychrom") ,]


#dxy_plot <- 
ggplot(dxy, aes(x=chrom, y=value, fill=chrom)) + guides(color = FALSE) +
  #geom_point(position=position_dodge(), stat="identity" ) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="Dxy") +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(. ~ sex) + guides(fill=FALSE) +
  scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#ffb14e', #orange
    '#1B75BB', #purple-y
    '#8B4BD8'
  ))

  annotate(geom="text", x = "A", y = 0.014, label = "a", parse = TRUE, size=10) +
  annotate(geom="text", x = "X", y = 0.006, label = "b", parse = TRUE, size=10) +
  annotate(geom="text", x = "Y", y = 0.007, label = "b", parse = TRUE, size=10) 
  
  t.test(
    all_data_var$value[all_data_var$chrom == "A" & all_data_var$var == "dxy" ] ,  
    all_data_var$value[all_data_var$chrom == "Xm" & all_data_var$var == "dxy" ])
  
  t.test(
    all_data_var$value[all_data_var$chrom == "A" & all_data_var$var == "dxy" ] ,  
    all_data_var$value[all_data_var$chrom == "Ychrom" & all_data_var$var == "dxy" ])
  
  t.test(
    all_data_var$value[all_data_var$chrom == "Ychrom" & all_data_var$var == "dxy" ] ,  
    all_data_var$value[all_data_var$chrom == "Xm" & all_data_var$var == "dxy" ])
  
  
multiplot(fst_plot,dxy_plot)


####kxy####
kxydxy <- rbind(
all_data[all_data$var == "kxy" & all_data$pop == "R.hastatulus",],
all_data[ all_data$var == "dxy",]
)
kxydxy$outgroup[kxydxy$var == "dxy"] <- "R.hastatulus"
kxydxy$outgroup[kxydxy$outgroup == "rothschildianus"] <- "R.rothschildianus"
kxydxy$outgroup[kxydxy$outgroup == "bucephalophorus"] <- "R.bucephalophorus"

kxydxy$Chromosome = factor(kxydxy$chrom, levels=c("A","N","H","X","Y"))


kxydxy$Outgroup =factor(kxydxy$outgroup, c("R.hastatulus","R.rothschildianus", "R.bucephalophorus"))


#levels(kxydxy$outgroup) = c("R.rothschildianus", "R.hastatulus")


#div_plot <- 
ggplot(kxydxy, aes(x=Chromosome, y=value, fill=Chromosome
                      # , group=1
                       )) + #guides(fill = FALSE) +
  #geom_point(position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + guides(fill=FALSE) +
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Average Divergence") +
 # scale_x_discrete(limits=c("R.hastatulus","R.rothschildianus","R.bucephalophorus")) +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  #theme(axis.text.x = element_text(face = "italic")) +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(sex ~ Outgroup) +
  scale_x_discrete(limits=c("A","N","H","X","Y")) +
  scale_fill_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#ffb14e', #orange
    '#1B75BB', #purple-y
    '#8B4BD8'
  ))
t.test(
  all_data_var$value[all_data_var$chrom == "A" & all_data_var$var == "dxy" ] ,  
  all_data_var$value[all_data_var$chrom == "X" & all_data_var$var == "dxy" ])

t.test(
  all_data_var$value[all_data_var$chrom == "Y" & all_data_var$var == "dxy" ] ,  
  all_data_var$value[all_data_var$chrom == "X" & all_data_var$var == "dxy" ])


dxy_data_var <- all_data_var[ all_data_var$var == "dxy", ]  
dxy_anova <- aov(value ~ chrom, data=dxy_data_var)
summary(dxy_anova) 



####dnds####
evo_rate <- rbind(
  all_data[all_data$var == "kaks" & all_data$pop == "R.hastatulus",],
  all_data[all_data$var == "dnds" ,]
  
)

#levels(all_data_dxy$outgroup) = c("rothschildianus", "hastatulus","bucephalophorus")
evo_rate$outgroup[evo_rate$var == "dnds"] <- "R.hastatulus"
evo_rate$outgroup[evo_rate$outgroup == "rothschildianus"] <- "R.rothschildianus"
evo_rate$outgroup[evo_rate$outgroup == "bucephalophorus"] <- "R.bucephalophorus"

evo_rate$Outgroup =factor(evo_rate$outgroup, c("R.hastatulus","R.rothschildianus", "R.bucephalophorus"))


#rate_plot <-
ggplot(evo_rate, aes(x=chrom, y=value, color=chrom
                         # , group=1
)) + #guides(fill = FALSE) +
  geom_point(position=position_dodge(.9)) +
  #geom_bar(position=position_dodge(), stat="identity" ) +
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Divergence (dn/ds)") +
  #scale_x_discrete(limits=c("R.hastatulus","R.rothschildianus","R.bucephalophorus")) +
  #theme(axis.text.x = element_text(face = "italic")) +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  theme(strip.background =element_rect(fill="white")) +
  scale_x_discrete(limits=c("A","N","H","X","Y")) + guides(color=FALSE) +
  
  facet_grid(sex ~ Outgroup) +
  scale_color_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#ffb14e', #orange
    '#1B75BB', #purple-y
    '#8B4BD8'
  ))


####mk####

mk <- all_data[all_data$var == "mk" & all_data$pop == "R.hastatulus",]
mk$outgroup[mk$outgroup == "rothschildianus"] <- "R.rothschildianus"
mk$outgroup[mk$outgroup == "bucephalophorus"] <- "R.bucephalophorus"


#mk$Chromosome = factor(mk$chrom, levels=c('A','XY','N','H'))
#mk$Chromosome = factor(mk$chrom, levels=c("A","N","H","X","Y"))

mk$Outgroup =factor(mk$outgroup, c("R.rothschildianus", "R.bucephalophorus"))


#adapt_plot <-
ggplot(mk, aes(x=chrom, y=value, color=chrom
                     # , group=1
)) + #guides(fill = FALSE) +
  geom_abline(intercept = 1,slope=0,alpha=0.2) +
  geom_point(position=position_dodge(.9)) +
  #geom_bar(position=position_dodge(), stat="identity" ) +
  geom_abline(intercept = 1,slope=0,alpha=0.2) +
  #stat_summary(fun.y=sum, geom="line") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 18) + labs(x="",y="Rate of Adaptation (MK)") +
  #scale_x_discrete(limits=c("R.rothschildianus","R.bucephalophorus")) +
  #theme(axis.text.x = element_text(face = "italic")) + coord_trans(y = "log2") +
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
  theme(strip.background =element_rect(fill="white")) +
  ylim(0.5,2) +
  scale_x_discrete(limits=c("A","N","H","X","Y")) + guides(color=FALSE) +
  facet_grid(sex ~ Outgroup) +
  scale_color_manual(values=c( 
    '#00ADEF', #Blue
    '#FFF100',  #yellow
    '#00A550', #green
    '#ffb14e', #orange
    '#1B75BB', #purple-y
    '#8B4BD8'
  ))

##test if group mean is diff from zero##
logmk_N <- data.frame(log(all_data_var$value[all_data_var$var == "mk" & all_data_var$pop == "R.hastatulus" & all_data_var$chrom == "N"]))
names(logmk_N) <- "logMK"
ggplot(logmk_N,aes(x=logMK)) + geom_histogram()
t.test(logmk_N)

####Qpi####

Q_data_pi <- data.frame( cbind(
    c("XY","XYY"),
    c(
    all_data_pi$value[all_data_pi$chrom == "XY" & all_data_pi$pop == "XY" & all_data_pi$outgroup == "rothschildianus"] / all_data_pi$value[all_data_pi$chrom == "A" & all_data_pi$pop == "XY" & all_data_pi$outgroup == "rothschildianus"],
    all_data_pi$value[all_data_pi$chrom == "XY" & all_data_pi$pop == "XYY" & all_data_pi$outgroup == "rothschildianus"] / all_data_pi$value[all_data_pi$chrom == "A" & all_data_pi$pop == "XYY" & all_data_pi$outgroup == "rothschildianus"]
    ) ) )
names(Q_data_pi) <- c("Pop","Qpi")

options(digits = 4)
Q_data_pi$Qpi<- type.convert(Q_data_pi$Qpi, na.strings = "NA", as.is = FALSE, dec = ".",
             numerals =  "no.loss")

title_qpi <- expression(paste("Q", ""[pi]))

ggplot(Q_data_pi, aes(x=Pop, y=Qpi)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_qpi) 

 ####Qfst####

Q_data_fst <-
  data.frame(cbind(
      "R.hastatulus",
      all_data_fst$value[all_data_fst$chrom == "XY" & all_data_fst$outgroup == "rothschildianus"] / all_data_fst$value[all_data_fst$chrom == "A" & all_data_fst$outgroup == "rothschildianus"]
    ))

names(Q_data_fst) <- c("Pop","Qfst")
Q_data_fst$Qfst<- type.convert(Q_data_fst$Qfst, na.strings = "NA", as.is = FALSE, dec = ".", numerals =  "no.loss")

title_qfst <- expression(paste("Q", ""[Fst]))

ggplot(Q_data_fst, aes(x=Pop, y=Qfst)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=title_qfst) 


####MS####


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
##fmaleX!=maleX##
maleX <- data.frame(stats_table(outgroup="rothschildianus",set="XYphased",chrom="X",pops=pop))
femaleX <- data.frame(stats_table(outgroup="rothschildianus",set="joshFem",chrom="X",pops=pop))

maleX[maleX$var == "pi",]
femaleX[femaleX$var == "pi",]

##synsite##
all_data[all_data$var == "sites" & all_data$cod == "syn" 
         & all_data$pop == "R.hastatulus" & all_data$chrom == "A", ]

#NE 

0.00629 / (4 * 7.5e-9) # Ne ancestral

(0.00629 / (4 * 7.5e-9) ) * 3.4 
(0.00629 / (4 * 7.5e-9) ) * 0.079

#Dxy = 2μt
0.269506300/(2*7.5e-9)
0.014862532/(2*7.5e-9)
(115000+35000)-0.014862532/(2*7.5e-9)

#Nem = ¼(fst - 1)  
(1/0.1165208 - 1)/4

#silene pollen expressed dnds
s_pollen_dnds <- 
data.frame(
  rbind(
c(0.192, 0.1689), 
c(0.142 ,0.1436)
), rownames = c("not pollen","pollen")
)
names(s_pollen_dnds) <- c("mean","se","tissue")

ggplot(s_pollen_dnds, aes(x=tissue, y=mean, color=tissue)) +
  guides(color = FALSE) +
  geom_point( ) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y="dnds") +
#theme(axis.text.x = element_text(angle = 20, hjust = 1))  +
#scale_x_discrete(limits=c("A","H","X","Y")) 

  scale_color_manual(values=c( 
    'red', 
    'blue'
    
  ))
