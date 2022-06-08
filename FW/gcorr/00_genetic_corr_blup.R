# Calculate genetic correlation between two traits
# Estimate genetic effect using rrBLUP 
# Contact: Hao Tong, tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "~/Nextcloud/05_fluxGxE/FW/"
setwd(dir)

library(rrBLUP)
options(digits = 15)

##########################################################################
## load data
# phenotypic data 
fw_blup <- read.table("FW-BLUP.csv",head=T,sep=",")[,2]
fw_blup <- matrix(fw_blup,ncol=1)
fcfw_blup <- read.table("FCFW-BLUP.csv",head=T,sep=",")[,2]
fcfw_blup <- matrix(fcfw_blup,ncol=1)

# genotypic data
accid <- read.table("../araid67_final_new.csv",head=F,sep=",")[,1]
myG <- read.table("../arasnp_num.csv",head=F,sep=",")
geno <- as.matrix(t(myG[,accid]))

# check accession names
#all(p2_id == snp_id[-1])

# normalization
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

fw_blup_norm <- apply(fw_blup,2,min_max_norm)
fcfw_blup_norm <- apply(fcfw_blup,2,min_max_norm)

pheno_all_norm <- cbind(fw_blup_norm,fcfw_blup_norm)

##########################################################################
### function rrBLUP model
rrBLUP <- function(geno,pheno,name){
  
  ucoef_all <- NULL
  for (n in 1:ncol(pheno)){
    Y <- as.numeric(pheno[,n])
    sol <- mixed.solve(Y,Z=geno,K=NULL,SE=F)
    ucoef <- as.matrix(sol$u)
    ucoef_all <- cbind(ucoef_all,ucoef)
  }
   
  write.table(ucoef_all,paste("gcorr/ucoef_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)

}

rrBLUP(geno,pheno_all_norm,"allblup")

### genetic correlation between traits
ucoef <- read.table("gcorr/ucoef_allblup.csv",head=F,sep=",")
corg <- cor(ucoef[,1],ucoef[,2],method="pearson",use="complete.obs")

rownames(corg) <- c("fw")
colnames(corg) <- c("fcfw")
write.table(corg,"gcorr/gcor_allblup.csv",sep=",",row.names=T,col.names=T)


