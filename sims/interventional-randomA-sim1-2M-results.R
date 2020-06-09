rm(list=ls())
libraries_check <- c("data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

setwd("sim1-2M/")
myfiles <- list.files()
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
sim_results_list <- NULL
for (ll in myfiles) {
  load(ll)
  seed <- as.integer(rev(strsplit(strsplit(ll,"-")[[1]],".Rdata"))[[1]])
  sim_results_list[[seed]] <- data.frame(do.call(rbind,lapply(simres, unlist)))
  rm(simres)
}

res <- rbindlist(sim_results_list)
setkey(res); rm(sim_results_list)

# number of sims & MC error
res[,.N,by=list(Um,corrM,contY,n,mc)]
# unique values of population effects
res[, round(unique(.SD),2),by=list(Um,corrM,contY),
    .SDcols=grep(pattern="pop",colnames(res))]

simset <- 4 # largest column index of sim setting

# direct effects
de <- res[, c(1:simset,grep("g0",names(res)),grep("t0",names(res))), with=FALSE]
de_summary <- de[,lapply(.SD, function(x) {
  paste0(formatC(mean(x),digits=2,format="f")," (",
         formatC(sd(x),digits=2,format="f"),")")
}),by=list(Um,corrM,contY,n)]
setkey(de_summary)
de_summary

# indirect effects
## interventional
ie <- res[, c(names(res)[1:4],
              grep("g1",names(res),value=TRUE), # joint indirect
              "ie.t1.pop", # indirect via M1 (true)
              "t1.intv", # indirect via M1 (est)
              "t2.intv.m2m1", # indirect via M1 (est; switched)
              "ie.t2.pop", # indirect via M2 (true)
              "t2.intv", # indirect via M2 (est)
              "t1.intv.m2m1", # indirect via M2 (est; switched)
              "ie.mu.pop", # indirect via mutual (true)
              grep("mu.intv",names(res),value=TRUE) # indirect via mutual (est)
              ),with=FALSE]
## natural
nie <- res[, c(1:simset,which(
  (grepl("natu",names(res)) & !grepl("[.]o[.]",names(res))) & 
    (grepl("t1",names(res)) | grepl("t2",names(res))))), 
  with=FALSE]
nie_types <- unique(unlist(
  lapply(strsplit(names(nie),"[.]t"),"[[",1)))[-(1:simset)]
nie_types
for (nn in 1:length(nie_types)) {
  nie_nn <- nie[,grep(paste0(nie_types[nn],".t"),names(nie)),with=FALSE]
  nie_nn[, paste0(nie_types[nn],".tsum") := rowSums(nie_nn)]
  setcolorder(nie_nn,c(4,1:3))
  if (grepl("m2m1",nie_types[nn])) {
    setcolorder(nie_nn,c(1,3,2,4)) # swap indirect effects for relabelled M
  }
  ie <- cbind(ie,nie_nn)
}
ie_summary <- ie[,lapply(.SD, function(x) {
  paste0(formatC(mean(x),digits=2,format="f")," (",
         formatC(sd(x),digits=2,format="f"),")")
}),by=list(Um,corrM,contY,n)]
setkey(ie_summary)

# no unobserved confounding of mediators, no mediator-mediator interaction
## all effects except for natural under incorrect direction should be unbiased
ie_summary[Um==0 & corrM==0 & contY==TRUE & n==500]
# unobserved confounding of mediators: natural effects should be biased
ie_summary[Um==1 & corrM==0 & contY==TRUE & n==500]
# correlated mediators differ by treatment: natural effects should be biased
ie_summary[Um==1 & corrM==1 & contY==TRUE & n==500]

# binary Y
ie_summary[contY==FALSE,c(1:simset,grep("pop",names(ie_summary))),with=FALSE]
library("xtable")
xtable(rbind(
  ie_summary[contY==FALSE, sort(which(
    !grepl("m2m1",names(ie_summary)) & !grepl("pop",names(ie_summary)))), 
    with=FALSE],
  ie_summary[contY==FALSE, c(1:simset,grep("m2m1",names(ie_summary))),with=FALSE],
  use.names=FALSE))

