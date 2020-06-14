rm(list=ls())
libraries_check <- c("data.table", "coxed", "xtable")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

source("hdmed-dataprep.R")

# load observed estimates from double selection ###############################
load(file="dataex-boot_ds-obs.Rdata")
any(is.na(res.ds)) # check no NAs

res.obs <- res.ds
rm(res.ds)

## order in Table of Huang and Pan (2016)
huangpan_table2 <- c("GO0009607","GO0051707","GO0016064",
                     "GO0019724","GO0045088","GO0002443",
                     "GO0002285","GO0002764","GO0006952")

# combine bootstrap samples from double selection #############################
if (!file.exists("interventional-dataex-2-ds-bootstrap-res.Rdata")) {
  subfolder <- "dataex-boot-ds/"
  myfiles <- list.files(subfolder)
  myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
  res_list <- NULL
  boot_ids <- rep(0,nrow(mydata1)); names(boot_ids) <- 1:nrow(mydata1)
  for (ll in myfiles) {
    # results
    loadable <- tryCatch(load(paste0(subfolder,ll)),
                         error=function(cond) return(NA))
    if(!is.na(loadable)) {
      if (exists("res.ds")) {
        res_list <- c(res_list,list(res.ds))
        for (i in boot_id) {
          boot_ids[i] <- boot_ids[i]+1
        }
        rm(boot_id,res,res.ds)
        cat(ll," boot ok \n")
      }
    }
  }
  plot(boot_ids)
  res_est <- do.call(cbind,res_list)
  dim(res_est)
  # check consecutive indices
  all(as.integer(lapply(strsplit(rownames(res_est),split="t"),"[",2))==
        1:nrow(res_est))
  all(as.integer(lapply(strsplit(names(res.obs),split="t"),"[",2))==
        1:nrow(res_est))
  save(res_est,file="interventional-dataex-2-ds-bootstrap-res.Rdata")
}

load(file="interventional-dataex-2-ds-bootstrap-res.Rdata")
# observed estimate (with bootstrap CIs)
res_est.ci <- cbind(res.obs,t(apply(res_est, 1, bca, conf.level=0.95)))
colnames(res_est.ci) <- c("est","l","u")

## replace with original names
rownames(res_est.ci) <- gene_names
head(res_est.ci)

## significant mediators in each gene set
res_est.ci.signif <- NULL
res_est.ci.signif.onetab <- NULL
for (gs in huangpan_table2) {
  gs_idx <- which(gene_names %in% as.character(all.geneset[[gs]][,"x"]))
  res_est.ci.gs <- data.frame(res_est.ci[gs_idx,])
  gs_signif <- which( (res_est.ci.gs[,"l"] > 0) | (res_est.ci.gs[,"u"] < 0) )
  res_est.ci.gs <- res_est.ci.gs[gs_signif,]
  
  res_est.ci.gs.dt <- data.table(cbind(
    "gene"=row.names(res_est.ci.gs),res_est.ci.gs,"gs"=1))
  setnames(res_est.ci.gs.dt,"gs",gs)
  setkey(res_est.ci.gs.dt)
  if (!is.null(res_est.ci.signif.onetab)) {
    res_est.ci.signif.onetab <- merge(res_est.ci.signif.onetab,
                                      res_est.ci.gs.dt,all=TRUE)
  } else {
    res_est.ci.signif.onetab <- res_est.ci.gs.dt
  }
  setkey(res_est.ci.signif.onetab)
}

res_est.ci.signif.onetab[is.na(res_est.ci.signif.onetab)] <- 0
print(xtable(res_est.ci.signif.onetab[order(res_est.ci.signif.onetab$l),],
             digits=c(0,0,2,2,2,rep(0,9))),include.rownames=FALSE)

# coefficients in outcome and mediator mean models ============================
coef.ds <- do.call(rbind,lapply(1:length(res), function(s) 
  c(res[[s]][["Y"]][paste0("M",s)],res[[s]][["M"]]["A"])))
row.names(coef.ds) <- gene_names
pdf("plot-dataex-iesets.pdf",width=9,height=9)
par(mfrow=c(3,3))
for (gs in huangpan_table2) {
  # genes in each set
  gs_idx <- which(gene_names %in% as.character(all.geneset[[gs]][,"x"]))
  maxColorValue <- 20
  palette <- colorRampPalette(colors=c("limegreen","red"),alpha=TRUE)(maxColorValue)
  # regression coefficients
  gs_regcoef <- coef.ds[gs_idx,]
  # indirect effects
  res_est.ci.gs <- data.frame(res_est.ci[gs_idx,])
  gs_ptsize <- pmin(exp(abs(res_est.ci.gs[,"est"]))*1.5,3)
  gs_signif <- which( (res_est.ci.gs[,"l"] > 0) | (res_est.ci.gs[,"u"] < 0) )
  
  plot(gs_regcoef,
       pch=20,
       col=palette[cut(gs_ptsize, maxColorValue)],
       cex=gs_ptsize,
       xlim=range(coef.ds[,1])*1.1,ylim=range(coef.ds[,2]),
       xlab="M --> Y",ylab="A --> M",
       main=gs)
  points(x=gs_regcoef[gs_signif,"M1"],y=gs_regcoef[gs_signif,"A"],
         cex=gs_ptsize)
  abline(h=0,col="grey50",lwd=.5)
  abline(v=0,col="grey50",lwd=.5)
}
dev.off()
