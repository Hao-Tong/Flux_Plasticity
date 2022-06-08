# Genome Selection using rrBLUP 
# Prediction accuracy and coincidence from cross-validations
# Contact: tong@mpimp-golm.mpg.de

# module load /apps/Modules/devel/R-4.1.2
# nohup Rscript 00_rrBLUP_new_LC.R > 00_rrBLUP_new_LC.out 2>&1 &

##########################################################################
## add path and packages
## module load /apps/Modules/devel/R-4.1.0
dir <- ("/winmounts/tong/winhome/mpidir/8.AlleleFlux/00_GxE/FC_LNLC/")
#dir <- ("H:/mpidir/8.AlleleFlux/00_GxE/FC_LNLC/")
#dir <- ("~/Nextcloud/05_fluxGxE/FC_LNLC/")
setwd(dir)
rm(list=ls())

library(rrBLUP)
options(digits = 15)

##########################################################################
## load data
# phenotypic data 
datafcln <- read.table("../FC/fcln_data.csv",head=F,sep=",")
datafclc <- read.table("../FC/fclc_data.csv",head=F,sep=",")
datafcln <- t(datafcln)
datafclc <- t(datafclc)

dataopt <- read.table("../opt/AraCore_allacc_opt.csv",sep=",",header=F)
dataln <- read.table("../LN/AraCore_allacc_LN.csv",sep=",",header=F)
datalc <- read.table("../LC/AraCore_allacc_LC.csv",sep=",",header=F)
dataopt <- t(dataopt)
dataln <- t(dataln)
datalc <- t(datalc)

# genotypic data
accid <- read.table("../araid67_final_new.csv",head=F,sep=",")[,1]
myG <- read.table("../arasnp_num.csv",head=F,sep=",")
geno <- as.matrix(t(myG[,accid]))

# cross-validation fold id
idsall <- read.table("fold/foldid.csv",sep=",",header=F)[,1:50]
rr <- 50 #replicate number
f <- 3 #fold number

##########################################################################
### function rrBLUP model
rrBLUP <- function(xx,name,pheno){
  
  for (n in 328:ncol(pheno)){
    
    corallall <- NULL
    ypred_all <- matrix(NA,nrow(pheno),rr)
    print(paste("Phenotype ",n," Start!",sep=""))
    Y <- as.numeric(pheno[,n])
    
    for (r in 1:rr){
      
      print(paste("Replicate ",r," Phenotype ",n," Done!",sep=""))
      ids <- idsall[,r]
      
      for (p in 1:f){
        trset <- which(ids!=p)
        #teset <- which(ids==p)
        xxtr <- xx[trset,]
        #xxte <- xx[teset,]
        ytr <- Y[trset]
        
        ff <- 3 #folds number
        nn <- length(trset) #sample number
        ny <- floor(nn/ff)
        id <- rep(c(1:ff),ny)
        if (nn > length(id)){
          id <- c(id,c(1:(nn-length(id))))
        }
        rrr <- 5 #replicate number
        foldall <- NULL
        for (r in 1:rrr){
          idss <- sample(id)
          foldall <- cbind(foldall,idss)
        }
        
        corall1 <- NULL
        for (zz in 1:rrr){
          ids1 <- foldall[,zz]
          
          for (pp in 1:ff){
            trset1 <- which(ids1!=pp)
            teset1 <- which(ids1==pp)
            xxtr1 <- xxtr[trset1,]
            xxte1 <- xxtr[teset1,]
            ytr1 <- ytr[trset1]
            yte1 <- ytr[teset1]
            
            sol <- mixed.solve(ytr1,Z=xxtr1,K=NULL,SE=F)
            ucoef <- as.matrix(sol$u)
            ypred1 <- rep(sol$beta,nrow(xxte1))+as.numeric(xxte1 %*% ucoef)
            optii <- dataopt[trset[teset1],n]
            cor1 <- cor(yte1*optii,ypred1*optii,method="pearson",use="complete.obs")
            corall1 <- c(corall1,cor1)
          }
        }
        
        corallall <- cbind(corallall,corall1)
    }
    
    }
    
    corallall <- rbind(corallall,colMeans(corallall,na.rm=T))
    write.table(corallall,paste("rrBLUP_within/pearson_",name,"_trait_",n,".csv",sep=""),sep=",",row.names=F,col.names=F)
}
}

rrBLUP(geno,"LC_within",datafclc)

# xx <- geno
# name <- "LN_within"
# pheno <- datafcln

