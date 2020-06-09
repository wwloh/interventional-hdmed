load("hdmed/biom12421-sup-0002-suppdata-code/Dataset_for_Biometrics.RData")
# id to preserve initial ordering of data
mydata1 <- data.table(cbind("id"=1:nrow(mydata1), mydata1))
setkey(mydata1)
nrow(mydata1)

# exposure for analysis ########################################################
# dichotomize exposure
mydata1[, A := as.integer(hsa.miR.223 > median(hsa.miR.223))]
setcolorder(mydata1,c(1,ncol(mydata1),2:(ncol(mydata1)-1)))
setkey(mydata1)

# outcome for analysis ########################################################
# compare suvival at 3 months (y) and death indicator
table("y"=mydata1$y,"death"=mydata1$death)

# loss to follow up (y = NA)
mydata1[is.na(y),list(y,fut,death)]
# compare between the treatment groups
mydata1[is.na(y),table(is.na(y),A)]
# recode as 0
mydata1[is.na(y),y := 0]

# compare suvival at 3 months (y) and death indicator, vs follow-up time
by(data=mydata1,INDICES=mydata1$y,FUN=function(x) summary(x$fut))
by(data=mydata1,INDICES=mydata1$death,FUN=function(x) summary(x$fut))

# mediators ###################################################################
#### check for NAs
colnames(mydata1)[1:13]
mydata1_genes <- mydata1[,11:ncol(mydata1)]
dim(mydata1_genes)
genes.noNA <- which(apply(mydata1_genes,2,function(x) all(!is.na(x))))
mydata1_genes <- mydata1_genes[,genes.noNA,with=FALSE]
dim(mydata1_genes)

gene_names <- colnames(mydata1_genes)
colnames(mydata1_genes) <- paste0("M",1:length(gene_names))

Lnames <- c("age","gender","race")
# propensity score model
fitA.form <- as.formula(paste0("A~",paste(Lnames,collapse="*")))

if (FALSE) {
# check for convergence and violations of positivity assumption
fitA <- glm(fitA.form, family = binomial("logit"), data = mydata1)
fitA.OK <- fitA$converged & !fitA$boundary &
  min(min(fitA$fitted.values),1-max(fitA$fitted.values))>.Machine$double.eps*1e1

pA1hat <- predict(fitA,type="response") # predicted prob. of treatment
pA0hat <- 1-pA1hat # predicted prob. of control
wt_A1 <- 1/pA1hat # inverse prob. of treatment weight
wt_A0 <- 1/pA0hat # inverse prob. of control weight
ind_A1 <- which(mydata1$A==1) # indices for individuals in treatment group
ind_A0 <- which(mydata1$A==0) # indices for individuals in control group
wt_A1.raw <- wt_A1[ind_A1]
wt_A1 <- wt_A1.raw/sum(wt_A1.raw)
wt_A0.raw <- wt_A0[ind_A0]
wt_A0 <- wt_A0.raw/sum(wt_A0.raw)
wt_df <- data.frame(rbind(cbind("A"=1,"wt"=wt_A1),
                          cbind("A"=0,wt_A0)))
wt_df[,1] <- as.factor(wt_df[,1])
library("ggplot2")
p <- ggplot(wt_df, aes(x=wt,colour=A,linetype=A)) + 
  geom_freqpoly(bins = 50, size=1.1)+
  theme_classic()+
  ggtitle("Sampling probabilities in each treatment group A")+
  xlab("Sampling probabilities")+
  theme(legend.position="top")+
  geom_vline(xintercept = 1/length(ind_A1), linetype=3)
ggsave(filename="plot-dataex-weights.pdf",plot=p,device="pdf",
       width=6,height=4,units="in")
}