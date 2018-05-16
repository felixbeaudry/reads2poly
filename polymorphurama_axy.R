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


####IMPORT####

y <-
  fread('goodsex_summarystats_Ychrom.txt')

x <-
  fread('goodsex_summarystats_Xchrom.txt')
x_inter <-
  fread('goodsex_interpop_Xchrom.txt')

a_raw <- 
  fread('goodsites_summarystats_may12.txt')
  #fread('summarystats_A_josh.txt')

locInfo <- fread('synteny_10322_hastatulus_transcripts.txt')
locInfo$buckwheat_position <- locInfo$buckwheat_position/1000000 

####Y####


ypol <- melt(y,id.vars = "locus")
ypols <- separate(ypol, variable, c("pop","var","cod"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
ypols_comp <- ypols[complete.cases(ypols), ]
ypols_comp$value <- as.numeric(ypols_comp$value)

ypols_comp$pop[ypols_comp$pop == "pop0"] <- "y"
ypols_comp$pop[ypols_comp$pop == "pop1"] <- "yTX"
ypols_comp$pop[ypols_comp$pop == "pop2"] <- "yNC"

ytgc <- summarySE(ypols_comp, measurevar="value", groupvars=c("var","pop","cod"))
yadf <- data.frame(ytgc)

ypi <- yadf[yadf$var == 'pi' & yadf$cod == "syn" & yadf$pop != "y",]
ytaj <- yadf[yadf$var == 'TajD',]
ytheta_syn_NC <- yadf[yadf$var == 'theta' & yadf$pop == "yNC" & yadf$cod == "syn",]
ytheta_syn_TX <- yadf[yadf$var == 'theta' & yadf$pop == "yTX" & yadf$cod == "syn",]
ytheta <- yadf[yadf$var == 'theta' & yadf$cod == "syn" & yadf$pop != "y",]

yfst <- yadf[yadf$var == 'fst'& yadf$cod == "syn",]
ydxy <- yadf[yadf$var == 'Dxy'& yadf$cod == "syn",]
ymk <- yadf[yadf$var == 'MK',]
yalpha <- yadf[yadf$var == 'alpha',]

####Y2####
Y2 <-
  fread('summarystats_Y_ncfl_mar21.xls')

Y2pol <- melt(Y2,id.vars = "locus")
Y2pols <- separate(Y2pol, variable, c("pop","var","cod"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
Y2pols_comp <- Y2pols[complete.cases(Y2pols), ]
Y2pols_comp$value <- as.numeric(Y2pols_comp$value)

Y2pols_comp$pop[Y2pols_comp$pop == "pop0"] <- "Y2"
Y2pols_comp$pop[Y2pols_comp$pop == "pop1"] <- "YNC"
Y2pols_comp$pop[Y2pols_comp$pop == "pop2"] <- "YFL"


Y2tgc <- summarySE(Y2pols_comp, measurevar="value", groupvars=c("var","pop","cod"))
Y2adf <- data.frame(Y2tgc)

Y2pi <- Y2adf[Y2adf$var == 'pi' & Y2adf$cod == "syn" & Y2adf$pop != "Y2",]
Y2taj <- Y2adf[Y2adf$var == 'TajD',]
Y2theta_syn_NC <- Y2adf[Y2adf$var == 'theta' & Y2adf$pop == "YNC" & Y2adf$cod == "syn",]
Y2theta_syn_FL <- Y2adf[Y2adf$var == 'theta' & Y2adf$pop == "YFL" & Y2adf$cod == "syn",]
Y2theta <- Y2adf[Y2adf$var == 'theta' & Y2adf$cod == "syn" & Y2adf$pop != "Y2",]


####X####
xpol <- melt(x,id.vars = "locus")
xpols <- separate(xpol, variable, c("pop","var","cod"), sep = "_", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
xpols_comp <- xpols[complete.cases(xpols), ]
xpols_comp$value <- as.numeric(xpols_comp$value)

xpols_comp$pop[xpols_comp$pop == "pop0"] <- "X"
xpols_comp$pop[xpols_comp$pop == "pop1"] <- "XTX"
xpols_comp$pop[xpols_comp$pop == "pop2"] <- "XNC"

xtgc <- summarySE(xpols_comp, measurevar="value", groupvars=c("var","pop","cod"))
xadf <- data.frame(xtgc)

7.280047e-03

xpi <- xadf[xadf$var == 'pi' & xadf$cod == "syn" & xadf$pop != "X",]
xtaj <- xadf[xadf$var == 'TajD',]
xtheta_syn_NC <- xadf[xadf$var == 'theta' & xadf$pop == "XNC" & xadf$cod == "syn",]
xtheta_syn_TX <- xadf[xadf$var == 'theta' & xadf$pop == "XTX" & xadf$cod == "syn",]
xtheta <- xadf[xadf$var == 'theta' & xadf$cod == "syn" & xadf$pop != "X",]

xfst <- xadf[xadf$var == 'fst'& xadf$cod == "syn",]
xdxy <- xadf[xadf$var == 'Dxy'& xadf$cod == "syn",]
xmk <- xadf[xadf$var == 'MK',]
xalpha <- xadf[xadf$var == 'alpha',]

####A####

a_all <- separate(a_raw, locus, c("1","locus","2"), sep = "_", remove = TRUE,
                  convert = FALSE, extra = "merge", fill = "left")

joined <- data.frame(sqldf('select a_all.*, locInfo.* from a_all left join locInfo on a_all.locus = locInfo.hastatulus_transcript'))

a <- joined[joined$TXjAuto == "1",]
a <- a[,c(2,4:42)]

#calculating Ne from dadi output
a_sites <- a[!is.na(a$pop0_sites_syn),]
sum(a_sites$pop0_sites_syn)
2508.81/sum(a_sites$pop0_sites_syn)
2508.81/(sum(a_sites$pop0_sites_syn)*7.5e-9*4)
(2508.81/(sum(a_sites$pop0_sites_syn)*7.5e-9*4))*0.3
##

pol <- melt(a,id.vars = "locus")
pols <- separate(pol, variable, c("pop","var","cod"), sep = "_", remove = TRUE,
                 convert = FALSE, extra = "merge", fill = "left")

pols_comp <- pols[complete.cases(pols), ]
pols_comp$value <- as.numeric(pols_comp$value)

pols_comp$pop[pols_comp$pop == "pop0"] <- "All"
pols_comp$pop[pols_comp$pop == "pop1"] <- "TX"
pols_comp$pop[pols_comp$pop == "pop2"] <- "NC"

tgc <- summarySE(pols_comp, measurevar="value", groupvars=c("var","pop","cod"))
adf <- data.frame(tgc)

pi <- adf[adf$var == 'pi' & adf$cod == 'syn' & adf$pop != "All", ,]
taj <- adf[adf$var == 'TajD',]
theta <- adf[adf$var == 'theta' & adf$pop != "All" & adf$cod == "syn",]


theta_syn_NC <- adf[adf$var == 'theta' & adf$pop == "NC" & adf$cod == "syn",]
theta_syn_TX <- adf[adf$var == 'theta' & adf$pop == "TX" & adf$cod == "syn",]

theta_syn_NC$value / 4 * 6.5e9 
theta_syn_TX$value / 4 * 6.5e9 

fst <- adf[adf$var == 'fst' & adf$cod == 'syn',]
dxy <- adf[adf$var == 'Dxy' & adf$cod == 'syn',]

alpha <- adf[adf$var == 'alpha',]

pi_anc = (pi$value[pi$pop == 'NC'] + pi$value[pi$pop == 'TX']) /2

da = dxy$value - pi_anc

####PolyRatio####

XAtheta_syn_NC <- xtheta_syn_NC$value / theta_syn_NC$value
XAtheta_syn_TX <- xtheta_syn_TX$value / theta_syn_TX$value
XA <- data.frame(c(theta_syn_NC$value,theta_syn_TX$value),c(xtheta_syn_NC$value,xtheta_syn_TX$value),c(XAtheta_syn_NC,XAtheta_syn_TX),c("XYY","XY"))
XA <- rename(XA, c("c.XAtheta_syn_NC..XAtheta_syn_TX."="theta", "c..XYY....XY.."="Populations"))


yAtheta_syn_NC <- ytheta_syn_NC$value / theta_syn_NC$value
yAtheta_syn_TX <- ytheta_syn_TX$value / theta_syn_TX$value
yA <- data.frame(c(theta_syn_NC$value,theta_syn_TX$value),c(ytheta_syn_NC$value,ytheta_syn_TX$value),c(yAtheta_syn_NC,yAtheta_syn_TX),c("XYY","XY"))
yA <- rename(yA, c("c.yAtheta_syn_NC..yAtheta_syn_TX."="yatheta", "c..XYY....XY.."="Populations"))


ggplot()

####ViolinGraphs####

a_pi_syn_popd <- pols_comp[pols_comp$var == 'pi' & pols_comp$cod == 'syn' & pols_comp$pop != "All", ,]
x_pi_syn_popd <- xpols_comp[xpols_comp$var == 'pi' & xpols_comp$cod == 'syn' & xpols_comp$pop != "X", ,]
y_pi_syn_popd <- ypols_comp[ypols_comp$var == 'pi' & ypols_comp$cod == 'syn' & ypols_comp$pop != "y", ,]
y2_pi_syn_popd <- Y2pols_comp[Y2pols_comp$var == 'pi' & Y2pols_comp$cod == 'syn' & Y2pols_comp$pop != "Y2", ,]

pi_syn_popd <- rbind(a_pi_syn_popd, x_pi_syn_popd,y_pi_syn_popd,y2_pi_syn_popd)
head(pi_syn_popd)

write.csv(pi_syn_popd, file = "pi_syn_popd.csv")
pi_syn_popd_e <-
  fread('pi_syn_popd_e.csv')

ggplot(pi_syn_popd_e) + 
  geom_violin(aes(x=Chromosome,y=pi_syn,fill=Population)) #+
  geom_point(aes(x=Chromosome,y=pi_syn,fill=Population,alpha=0.5)) 

ggplot(pi_syn_popd_e, aes(x=Chromosome, y=pi_syn, fill=Population)) + 
    geom_violin(trim=FALSE)+
    # geom_boxplot(width=0.1, fill="white") +
theme_classic()

####BarGraphs####

theta_all <- rbind(ytheta, Y2theta, xtheta, theta)
theta_all["Location"] <- c("Y","Y","Ysub","Ysub","X","X","A","A")
theta_all["Population"] <- c("XYY","XY","FL","SC","XYY","XY","XYY","XY")

pi_all <- rbind(pi,ypi, xpi)
pi_all["Location"] <- c("a","a","y","y","x","x")
pi_all["Population"] <- c("XYY","XY","XYY","XY","XYY","XY")

fst_all <- rbind(fst, xfst,  yfst)
fst_all["Location"] <- c("A","X", "Y")

dxy_all <- rbind(dxy, xdxy, ydxy)
dxy_all["Location"] <- c("A","X","Y")

alpha_all <- rbind(yalpha, alpha )
alpha_all["Location"] <- c("y","y","A","A")
alpha_all["Population"] <- c("XYY","XY","XYY","XY")

title3 <- expression(paste(theta, ""[syn]))

#thetaplot <- 
  ggplot(theta_all, aes(x=Location, y=value, fill=Population)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  theme_bw()  + theme_bw(base_size = 30) +
  labs(x = " ",y = title3)  # + scale_fill_manual(values=c('light blue','dark blue')) +
  
  annotate(geom="text", x = 0.77, y = 0.012, label = "a", parse = TRUE, size=15) +
  annotate(geom="text", x = 1.23, y = 0.009, label = "b", parse = TRUE, size=15) +
  annotate(geom="text", x = 1.77, y = 0.007, label = "c", parse = TRUE, size=15) +
  annotate(geom="text", x = 2.23, y = 0.0045, label = "d", parse = TRUE, size=15) 

  theta_ys <- rbind(ytheta, Y2theta)
  theta_ys["Location"] <- c("Y","Y","Ysub","Ysub")
  theta_ys["Population"] <- c("XYY","XY","FL","SC")
  
  theta_ys <- theta_ys[-c(1),]
  
  ggplot(theta_ys, aes(x=Location, y=value, fill=Population)) + #guides(fill = FALSE) +
    geom_bar(position=position_dodge(), stat="identity" ) +
    geom_errorbar(aes(ymin=value-se, ymax=value+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
    theme_bw()  + theme_bw(base_size = 30) 
  
  
  
  
  pi_ys <- rbind(ypi, Y2pi )
  pi_ys["Location"] <- c("Y","Y","Ysub","Ysub")
  pi_ys["Population"] <- c("XYY","XY","FL","SC")
  
  pi_ys <- pi_ys[-c(1),]
  
  titlepi <- expression(paste(pi, ""[syn]))
  
  ggplot(pi_all, aes(x=Location, y=value, fill=Population)) + #guides(fill = FALSE) +
    geom_bar(position=position_dodge(), stat="identity" ) +
    geom_errorbar(aes(ymin=value-se, ymax=value+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
    theme_bw()  + theme_bw(base_size = 30) + labs(x = "", y=titlepi)
  
t.test(
  Y2pols_comp$value[Y2pols_comp$pop == "YFL" & Y2pols_comp$var == "pi" & Y2pols_comp$cod == "syn"],  
  Y2pols_comp$value[Y2pols_comp$pop == "YNC" & Y2pols_comp$var == "pi" & Y2pols_comp$cod == "syn"])
#p-value = 0.1213


#combine y2 with pollen bias



Y2_split <- separate(Y2, locus, c("1","locus","2"), sep = "_", remove = TRUE,
                             convert = FALSE, extra = "merge", fill = "left")


y2_pisyn_pollen <- data.frame(sqldf('select gFPKM_pf_split.locus, gFPKM_pf_split.P2F, gFPKM_pf_split.V2,
                                      Y2_split.pop1_pi_syn,  Y2_split.pop2_pi_syn
  
                                       from gFPKM_pf_split 
                                      left join Y2_split on  Y2_split.locus = gFPKM_pf_split.locus'))


y2_pisyn_pollen_comp <- y2_pisyn_pollen[complete.cases(y2_pisyn_pollen), ]

y2_pisyn_pollen_melt <- 
  melt(y2_pisyn_pollen_comp,id.vars = c("locus","V2","P2F"))


y2_pisyn_pollen_summary <- 
  data.frame(summarySE(y2_pisyn_pollen_melt, measurevar="value", groupvars=c("V2","variable")))



y2_pisyn_pollen_summary <- cbind(y2_pisyn_pollen_summary,pop=c("SC","FL","SC","FL","SC","FL"))


ggplot(y2_pisyn_pollen_summary, aes(x=V2, y=value, fill=pop)) + #guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  
  theme_bw()  + theme_bw(base_size = 30) + 
  labs(x = " ",y = titlepi)   + scale_fill_manual(values=c('light blue','dark blue')) 


t.test(
  y2_pisyn_pollen_melt$value[y2_pisyn_pollen_melt$V2 == "Pollen" & y2_pisyn_pollen_melt$variable == "pop1_pi_syn"] ,  
  y2_pisyn_pollen_melt$value[y2_pisyn_pollen_melt$V2 == "Pollen" & y2_pisyn_pollen_melt$variable == "pop2_pi_syn"])

fit <- aov(value ~ V2 + variable, data=y2_pisyn_pollen_melt)

summary(fit) # display Type I ANOVA table


#pi-n/pi-s
fl_pi_rs <- data.frame(cbind(Y2$pop1_pi_rep,Y2$pop1_pi_syn),3)
fl_pi_rs$pnps <- fl_pi_rs$X1/fl_pi_rs$X2
fl_pi_rs$pnps[is.infinite(fl_pi_rs$pnps)] <- (fl_pi_rs$X1[is.infinite(fl_pi_rs$pnps)]) + 1
fl_pi_rs$pnps[is.nan(fl_pi_rs$pnps)] <- 1
fl_pi_rs$pnps[fl_pi_rs$pnps == 0 ] <- 1 - (fl_pi_rs$X2[fl_pi_rs$pnps == 0] )

sc_pi_rs <-  data.frame(cbind(as.numeric(Y2$pop2_pi_rep),as.numeric(Y2$pop2_pi_syn),2))
sc_pi_rs$pnps <- sc_pi_rs$X1/sc_pi_rs$X2
sc_pi_rs$pnps[is.infinite(sc_pi_rs$pnps)] <- (sc_pi_rs$X1[is.infinite(sc_pi_rs$pnps)]) + 1
sc_pi_rs$pnps[is.nan(sc_pi_rs$pnps)] <- 1
sc_pi_rs$pnps[sc_pi_rs$pnps == 0 ] <- 1 - (sc_pi_rs$X2[sc_pi_rs$pnps == 0])

y2_pnps <- rbind(fl_pi_rs,sc_pi_rs)

y2_pnps_sum <- 
  summarySE(y2_pnps, measurevar="pnps", groupvars=c("X3"))

y2_pnps_sum$X3[y2_pnps_sum$X3 == 2] <- "SC"
y2_pnps_sum$X3[y2_pnps_sum$X3 == 3] <- "FL"

titlepin <- expression(paste(pi, ""[rep],"/",pi, ""[syn]))

ggplot(y2_pnps_sum, aes(x=X3, y=pnps, fill=X3)) + guides(fill = FALSE) +
  geom_bar(position=position_dodge(), stat="identity" ) +
  geom_errorbar(aes(ymin=pnps-se, ymax=pnps+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  
  theme_bw()  + theme_bw(base_size = 30) + ylim(0.8,1.2) +
  labs(x = " ",y = titlepin)   + scale_fill_manual(values=c('light blue','dark blue')) 



t.test(sc_pi_rs$pnps,fl_pi_rs$pnps)
#0.3157


####SlidingGraphs####

#y
y_sep <- separate(y, locus, c("1","locus","2"), sep = "_", remove = TRUE,convert = FALSE, extra = "merge", fill = "left")
y_joined <- data.frame(sqldf('select y_sep.*, locInfo.* from y_sep left join locInfo on y_sep.locus = locInfo.hastatulus_transcript'))
#y_pisyn_pos <- y_joined[,c(19,30,46,47)]
#y_pisyn_pos <- rename(y_pisyn_pos, c("pop1_pi_syn"="NC", "pop2_pi_syn"="TX"))
y_pisyn_pos <- y_joined[,c(30,46,47)]
y_pisyn_pos <- rename(y_pisyn_pos, c( "pop2_pi_syn"="yTX"))
y_pisyn_pos <- y_pisyn_pos[y_pisyn_pos$buckwheat_chromosome < 9, ]
y_pisyn_melt <- melt(y_pisyn_pos,id.vars = c("buckwheat_chromosome","buckwheat_position"))
y_pisyn_melt <-  y_pisyn_melt[complete.cases(y_pisyn_melt), ]


y2_sep <- separate(Y2, locus, c("1","locus","2"), sep = "_", remove = TRUE,convert = FALSE, extra = "merge", fill = "left")
y2_joined <- data.frame(sqldf('select y2_sep.*, locInfo.* from y2_sep left join locInfo on y2_sep.locus = locInfo.hastatulus_transcript'))
y2_pisyn_pos <- y2_joined[,c(19,30,46,47)]
y2_pisyn_pos <- rename(y2_pisyn_pos, c("pop1_pi_syn"="ySC", "pop2_pi_syn"="yFL"))

y2_pisyn_pos <- y2_pisyn_pos[y2_pisyn_pos$buckwheat_chromosome < 9, ]
y2_pisyn_melt <- melt(y2_pisyn_pos,id.vars = c("buckwheat_chromosome","buckwheat_position"))
y2_pisyn_melt <-  y2_pisyn_melt[complete.cases(y2_pisyn_melt), ]

ys_pisyn <- rbind(y_pisyn_melt, y2_pisyn_melt)


#x
x_sep <- separate(x, locus, c("1","locus","2"), sep = "_", remove = TRUE,convert = FALSE, extra = "merge", fill = "left")
x_joined <- data.frame(sqldf('select x_sep.*, locInfo.* from x_sep left join locInfo on x_sep.locus = locInfo.hastatulus_transcript'))
x_pisyn_pos <- x_joined[,c(19,30,46,47)]
x_pisyn_pos <- rename(x_pisyn_pos, c("pop1_pi_syn"="xNC", "pop2_pi_syn"="xTX"))
x_pisyn_pos <- x_pisyn_pos[x_pisyn_pos$buckwheat_chromosome < 9, ]
x_pisyn_melt <- melt(x_pisyn_pos,id.vars = c("buckwheat_chromosome","buckwheat_position"))
x_pisyn_melt <-  x_pisyn_melt[complete.cases(x_pisyn_melt), ]


#a

a_pisyn_pos <- a[,c(19,30,68,69)]
a_pisyn_pos <- rename(a_pisyn_pos, c("pop1_pi_syn"="aNC", "pop2_pi_syn"="aTX"))
a_pisyn_pos <- a_pisyn_pos[a_pisyn_pos$buckwheat_chromosome < 9, ]
a_pisyn_melt <- melt(a_pisyn_pos,id.vars = c("buckwheat_chromosome","buckwheat_position"))
a_pisyn_melt <-  a_pisyn_melt[complete.cases(a_pisyn_melt), ]

title <- expression(paste( pi, ""[syn]))

a_pisyn_melt7 <- a_pisyn_melt[a_pisyn_melt$buckwheat_chromosome == 7,]

a_plot <-
  ggplot(a_pisyn_melt7,aes(x=buckwheat_position,y=value)) + 
  geom_point(aes(color=variable) ) +
  geom_smooth(aes(color=variable,alpha=0.7),se = FALSE,method="loess",span=1.5,size=2) + 
  facet_grid(. ~ buckwheat_chromosome, scales = "free") +
  labs(x = "",y = title,title="Autosomal") + #ylim(-0,1) +
  theme_bw(base_size = 18)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(size=FALSE,alpha=FALSE) #+ scale_color_manual(values=c( "#009E73",    "#E69F00","#56B4E9"))

x_pisyn_melt7 <- x_pisyn_melt[x_pisyn_melt$buckwheat_chromosome == 7,]

x_plot <-
  ggplot(x_pisyn_melt7,aes(x=buckwheat_position,y=value)) + 
  geom_smooth(aes(color=variable,alpha=0.7),se = FALSE,method="loess",span=1.5,size=2) + 
  geom_point(aes(color=variable) ) +
  facet_grid(. ~ buckwheat_chromosome, scales = "free") +
  labs(x = "",y = title,title="X") + #ylim(-0,1) +
  theme_bw(base_size = 18)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(size=FALSE,alpha=FALSE) #+ scale_color_manual(values=c( "#009E73",    "#E69F00","#56B4E9"))

#y_plot <-
  ggplot(ys_pisyn,aes(x=buckwheat_position,y=value)) + 
    geom_point(aes(color=variable,alpha=0.5) ) +
  geom_smooth(aes(color=variable,alpha=0.7),se = FALSE,method="loess",span=1.5,size=2) + 
  
  facet_grid(. ~ buckwheat_chromosome, scales = "free") +
  labs(x = "Position (MB)",y = title,title="Y") + #ylim(-0,1) +
  theme_bw(base_size = 18)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(size=FALSE,alpha=FALSE) #+ scale_color_manual(values=c( "#009E73",    "#E69F00","#56B4E9"))


multiplot(a_plot,x_plot)


#boxplot#
AXY_pisyn <- rbind(y_pisyn_melt[,c(1,3,4)], y2_pisyn_melt[,c(1,3,4)],x_pisyn_melt[,c(1,3,4)],a_pisyn_melt[,c(1,3,4)])

AXY_pisyn_stats <- summarySE(AXY_pisyn, measurevar="value", groupvars=c("buckwheat_chromosome","variable"))

write.csv(AXY_pisyn_stats, file = "AXY_pisyn_stats.csv")
AXY_pisyn_stats_e <-
  fread('AXY_pisyn_stats_e.csv')

AXY_pisyn_stats_e$chrom[AXY_pisyn_stats_e$chrom == "a"] <- "A"
AXY_pisyn_stats_e$chrom[AXY_pisyn_stats_e$chrom == "x"] <- "X"
AXY_pisyn_stats_e$chrom[AXY_pisyn_stats_e$chrom == "y"] <- "Y"

ggplot(AXY_pisyn_stats_e,aes(x=chrom,y=value,fill=pop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  facet_grid(. ~ buckwheat_chromosome, scales = "free") + labs(x = "Chromosome",y = title) +
  theme_bw(base_size = 18)  +  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
               width=.2,                    # Width of the error bars
             position=position_dodge(.9)) +
  scale_fill_manual(values=c( "#ffd700",
                               "#ffb14e",
                               "#fa8775",
                               "#ea5f94"))

AXY_pisyn_pos <- rbind(y_pisyn_melt, y2_pisyn_melt,x_pisyn_melt,a_pisyn_melt)

AXY_pisyn_pos_7 <- AXY_pisyn_pos[AXY_pisyn_pos$buckwheat_chromosome == 7, ]

max7 <- as.numeric(max(AXY_pisyn_pos_7$buckwheat_position)/3)

AXY_pisyn_pos_7$Position[AXY_pisyn_pos_7$buckwheat_position <= max7 ] <- "<17MB"
AXY_pisyn_pos_7$Position[AXY_pisyn_pos_7$buckwheat_position > max7 & AXY_pisyn_pos_7$buckwheat_position <= 2*max7 ] <- ">17MB & <34MB"
AXY_pisyn_pos_7$Position[AXY_pisyn_pos_7$buckwheat_position > 2*max7] <- ">34MB"

AXY_pisyn_stats_7 <- summarySE(AXY_pisyn_pos_7[,-c(1)], measurevar="value", groupvars=c("Position","variable"))

ggplot(AXY_pisyn_stats_7,aes(x=variable,y=value,fill=Position)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  #facet_grid(. ~ buckwheat_chromosome, scales = "free") + 
  labs(x = "",y = title,title="LG7") +
  theme_bw(base_size = 18)  +  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) #+
  scale_fill_manual(values=c( "#ffd700",
                              "#ffb14e",
                              "#fa8775",
                              "#ea5f94"))

