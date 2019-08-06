suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(sqldf))

subsetSFS<-function(fz=NA,list=NA,listString="subset",numInds=NA,sumstats=NA,pop = NA,sampSubSize=0.8){
  foldSFS <- function(i=NA){
    fold <- rep.int(0, length(i))
    fold[1] <- i[1]
    for (site in c(1:(round((length(i)-1)/2)))){
      fold[site+1] <-  i[site+1] + i[length(i)+1-site]
    }
    if( ( (length(i) - 1) %% 2) != 0 ){
      remainder <- (length(i)  / 2) +1 
      fold[remainder] <- i[remainder]
    }
    return(fold)
  }
  findInvar <- function(fz=NA,sumstats=NA, pop = pop, site = "syn"){
   
    numInds = length(fz) 
    numColsTemp <- c(2:numInds)
    fz$sums <- rowSums(fz[,..numColsTemp], na.rm = FALSE, dims = 1)
    fz_sep <- separate(fz, V1, c("locus","invar"), sep = ".fas", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
    sqlString <- paste("select fz_sep.*, sumstats.pop",pop,"_sites_",site," AS site, sumstats.pop",pop,"_pi_rate AS rate from fz_sep left join sumstats on fz_sep.locus = sumstats.locus",sep="")
    fz_sql <- data.frame(sqldf(sqlString))
    
    fz_sql$invar <- round(fz_sql$site - fz_sql$sums)
    
    return(fz_sql)
  }
  
  numColsSyn <- c(3:(numInds+1))
  numColsRep <- c((numInds+4):(2*numInds+2))
  
  synfz <- findInvar(fz=allelefz[,c(1,..numColsSyn)],sumstats=sumstats, pop = pop, site = "syn")
  nsfz <- findInvar(fz=allelefz[,c(1,..numColsRep)],sumstats=sumstats, pop = pop, site = "rep")
  
  sAll <- synfz[synfz$locus %in% list,]
  
  sample_tot <- length(sAll$locus)
  sample_sub <- sample(1:sample_tot, sample_tot*sampSubSize, replace=FALSE)
  
  s <- sAll[sample_sub, ]
  nsAll <- nsfz[nsfz$locus %in% list,]
  ns <- nsAll[sample_sub, ]
  
  do.call(cat,list(c("1")))
  do.call(cat,list(c("\n")))
  do.call(cat,list(c(numInds-1)))
  do.call(cat,list(c("\n")))
  
  ns_fold <- foldSFS(i= colSums(ns[,c(2:(numInds+1))], na.rm = FALSE, dims = 1))
  do.call(cat,list(ns_fold))
  ns_fz <- ns_fold/sum(ns_fold)
  
  do.call(cat,list(c("\n")))
  
  s_fold <- foldSFS(i=colSums(s[,c(2:(numInds+1))], na.rm = FALSE, dims = 1))
  do.call(cat,list(s_fold))
  s_fz <- s_fold/sum(s_fold)
  
  p<- data.frame(t(rbind.data.frame( t( c(0:(numInds-1)) ),s_fz, ns_fz )  ))
  names(p) <- c("frequency","Synonymous","NonSyn")
  pmelt <- melt(p,id.vars = "frequency",verbose=FALSE)
  pmelt <- cbind(pmelt,subD = listString)
  return(pmelt)
}

#import
sumstats <- fread('na_rothschildianus_summarystats_pop_X.txt',header=TRUE)
allelefz <- fread('na_rothschildianus_frequencies_pop1_X.txt',header=FALSE)

sumstats_sep <- separate(sumstats, locus, c("locus","butt"), 
                         sep = ".fas", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")

args <- commandArgs(trailingOnly = TRUE)
list <- fread(args,header = FALSE)

#list <- fread('pollen.list',header = FALSE)
list <- list$V1

#SFS <-  subsetSFS(fz=allelefz,list=list,listString="Yes",numInds=14,sumstats=sumstats_sep,pop="1",sampSubSize=0.8)

##alpha
sumstats_out <- sumstats_sep[!is.na(sumstats_sep$pop1_k_syn) & !is.na(sumstats_sep$pop1_k_rep) & sumstats_sep$locus %in% list,]
list <- sumstats_sep$locus[!is.na(sumstats_sep$pop1_k_syn) & !is.na(sumstats_sep$pop1_k_rep) & sumstats_sep$locus %in% list]

do.call(cat,list(c("1",sum(sumstats_out$pop1_sites_rep),sum(sumstats_out$pop1_k_rep*sumstats_out$pop1_sites_rep))))
do.call(cat,list(c("\n")))
do.call(cat,list(c("0",sum(sumstats_out$pop1_sites_syn),sum(sumstats_out$pop1_k_syn*sumstats_out$pop1_sites_syn))))



