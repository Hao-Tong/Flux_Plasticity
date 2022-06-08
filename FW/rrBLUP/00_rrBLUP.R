# Genome Selection using rrBLUP 
# Prediction accuracy and coincidence from cross-validations
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
#dir <- ("/winmounts/tong/winhome/mpidir/8.AlleleFlux/00_GxE/FC")
dir <- ("D:/05_fluxGxE/FW/")
setwd(dir)
rm(list=ls())

library(rrBLUP)
options(digits = 15)

##########################################################################
## load data
# phenotypic data 
fwdata <- read.table("../araid67_final_new_env.csv",head=F,sep=",")[,2:4]
pheno_opt <- as.matrix(fwdata[,1])
pheno_ln <- as.matrix(fwdata[,2])
pheno_lc <- as.matrix(fwdata[,3])

pheno_fcln <- pheno_ln/pheno_opt
pheno_fclc <- pheno_lc/pheno_opt

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
  
  for (n in 1:ncol(pheno)){
      
    ypred_all <- matrix(NA,nrow(pheno),rr)
    print(paste("Phenotype ",n," Start!",sep=""))
    Y <- as.numeric(pheno[,n])
    
    for (r in 1:rr){
      
      print(paste("Replicate ",r," Phenotype ",n," Done!",sep=""))
      ids <- idsall[,r]
      
      for (p in 1:f){
        trset <- which(ids!=p)
        teset <- which(ids==p)
        xxtr <- xx[trset,]
        xxte <- xx[teset,]
        ytr <- Y[trset]
        sol <- mixed.solve(ytr,Z=xxtr,K=NULL,SE=F)
        ucoef <- as.matrix(sol$u)
        ypred <- rep(sol$beta,nrow(xxte))+as.numeric(xxte %*% ucoef)
        ypred_all[teset,r] <- ypred
      }
    }
    
    write.table(ypred_all,paste("rrBLUP/predict_",name,"_trait_",n,".csv",sep=""),sep=",",row.names=F,col.names=F)
    
  }
}

### function predictability of regression models
ability_regression <- function(namein,pheno,nameout){
  
  cor1all <- matrix(NA,rr*f,ncol(pheno))
  cor2all <- matrix(NA,rr*f,ncol(pheno))
  # match1all <- matrix(NA,rr*f,ncol(pheno))
  # match2all <- matrix(NA,rr*f,ncol(pheno))
  
  for (n in 1:ncol(pheno)){
    
    print(paste("Phenotype ",n," Start!",sep=""))
    Ypred <- read.table(paste("rrBLUP/predict_",namein,"_trait_",n,".csv",sep=""),head=F,sep=",")
    Y <- as.numeric(pheno[,n])
    
    for (r in 1:rr){
      
      ids <- idsall[,r]
      cor1p <- NULL
      cor2p <- NULL
      #match1p <- NULL
      #match2p <- NULL
      
      for (p in 1:f){
        teset <- which(ids==p)
        ypred <- Ypred[teset,r]
        yte <- Y[teset]
      
        # accuracy
        cor1 <- cor(yte,ypred,method="pearson",use="complete.obs")
        cor2 <- cor(yte,ypred,method="spearman",use="complete.obs")
        cor1p <- c(cor1p,cor1)
        cor2p <- c(cor2p,cor2)
        
        # # coincidence
        # ytetop <- order(yte,decreasing=T)[1:round(length(yte)*0.3)] #top 30% genotypes
        # ypredtop <- order(ypred,decreasing=T)[1:round(length(yte)*0.3)]
        # match1 <- sum(ytetop %in% ypredtop)/length(ytetop)
        # match1p <- c(match1p,match1)
        # 
        # # coincidence
        # ytetop <- order(yte,decreasing=T)[1:round(length(yte)*0.15)] #top 15% genotypes
        # ypredtop <- order(ypred,decreasing=T)[1:round(length(yte)*0.15)]
        # match2 <- sum(ytetop %in% ypredtop)/length(ytetop)
        # match2p <- c(match2p,match2)
      }
      
      cor1all[(f*r-(f-1)):(f*r),n] <- cor1p
      cor2all[(f*r-(f-1)):(f*r),n] <- cor2p
      # match1all[(f*r-(f-1)):(f*r),n] <- match1p
      # match2all[(f*r-(f-1)):(f*r),n] <- match2p
    }
    
  }
  
  cor1allm <- rbind(cor1all,colMeans(cor1all,na.rm=T))
  cor2allm <- rbind(cor2all,colMeans(cor2all,na.rm=T))
  # match1allm <- rbind(match1all,colMeans(match1all,na.rm=T))
  # match2allm <- rbind(match2all,colMeans(match2all,na.rm=T))
  
  # Output: prediction accuracy and coincidence
  write.table(cor1allm,paste("rrBLUP/pearson_",nameout,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(cor2allm,paste("rrBLUP/spearman_",nameout,".csv",sep=""),sep=",",row.names=F,col.names=F)
  # write.table(match1allm,paste("rrBLUP/coincidence30_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  # write.table(match2allm,paste("rrBLUP/coincidence15_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)

  results <- rbind(colMeans(cor1all,na.rm=T),colMeans(cor2all,na.rm=T),
                   apply(cor1all,2,sd,na.rm=T),apply(cor2all,2,sd,na.rm=T))
  results <- t(results)
  results <- cbind("FW",results)
  colnames(results) <- c("ID","Mean-Pearson","Mean-Spearman","SD-Pearson","SD-Spearman")
  write.table(results,paste("rrBLUP/summary_",nameout,".csv",sep=""),sep=",",row.names=F,col.names=T)

}

##########################################################################
## FW in optimal
rrBLUP(geno,"opt",pheno_opt)
ability_regression("opt",pheno_opt,"opt2opt")
ability_regression("opt",pheno_ln,"opt2ln")
ability_regression("opt",pheno_lc,"opt2lc")
ability_regression("opt",pheno_fcln,"opt2fcln")
ability_regression("opt",pheno_fclc,"opt2fclc")

## FW in FC LN
rrBLUP(geno,"ln",pheno_ln)
ability_regression("ln",pheno_ln,"ln2ln")
ability_regression("ln",pheno_fcln,"ln2fcln")

rrBLUP(geno,"fcln",pheno_fcln)
ability_regression("fcln",pheno_fcln,"fcln2fcln")
ability_regression("fcln",pheno_opt,"fcln2opt")
ability_regression("fcln",pheno_ln,"fcln2ln")

## FW in FC LC
rrBLUP(geno,"lc",pheno_lc)
ability_regression("lc",pheno_lc,"lc2lc")
ability_regression("lc",pheno_fclc,"lc2fclc")

rrBLUP(geno,"fclc",pheno_fclc)
ability_regression("fclc",pheno_fclc,"fclc2fclc")
ability_regression("fclc",pheno_opt,"fclc2opt")
ability_regression("fclc",pheno_lc,"fclc2lc")

# FW in FC LN * raw optimal
Ypred_opt <- pheno_opt
Ypred_fcln <- read.table("rrBLUP/predict_fcln_trait_1.csv",head=F,sep=",")

Ypred_fcln.rawopt <- matrix(0,nrow(Ypred_fcln),ncol(Ypred_fcln))
for (m in 1:ncol(Ypred_fcln)){
  Ypred_fcln.rawopt[,m] <- as.numeric(Ypred_opt) * as.numeric(Ypred_fcln[,m])
}
write.table(Ypred_fcln.rawopt,"rrBLUP/predict_fcln.rawopt_trait_1.csv",sep=",",row.names=F,col.names=F)

ability_regression("fcln.rawopt",pheno_ln,"fcln.rawopt2ln")

# FW in FC LC * raw optimal
Ypred_opt <- pheno_opt
Ypred_fclc <- read.table("rrBLUP/predict_fclc_trait_1.csv",head=F,sep=",")

Ypred_fclc.rawopt <- matrix(0,nrow(Ypred_fclc),ncol(Ypred_fclc))
for (m in 1:ncol(Ypred_fclc)){
  Ypred_fclc.rawopt[,m] <- as.numeric(Ypred_opt) * as.numeric(Ypred_fclc[,m])
}
write.table(Ypred_fclc.rawopt,"rrBLUP/predict_fclc.rawopt_trait_1.csv",sep=",",row.names=F,col.names=F)

ability_regression("fclc.rawopt",pheno_lc,"fclc.rawopt2lc")

# FW in FC LN * CV optimal
Ypred_opt <- read.table("rrBLUP/predict_opt_trait_1.csv",head=F,sep=",")
Ypred_fcln <- read.table("rrBLUP/predict_fcln_trait_1.csv",head=F,sep=",")
Ypred_fcln.opt <- Ypred_opt * Ypred_fcln
write.table(Ypred_fcln.opt,"rrBLUP/predict_fcln.opt_trait_1.csv",sep=",",row.names=F,col.names=F)
ability_regression("fcln.opt",pheno_ln,"fcln.opt2ln")

# FW in FC LC * CV optimal
Ypred_opt <- read.table("rrBLUP/predict_opt_trait_1.csv",head=F,sep=",")
Ypred_fclc <- read.table("rrBLUP/predict_fclc_trait_1.csv",head=F,sep=",")
Ypred_fclc.opt <- Ypred_opt * Ypred_fclc
write.table(Ypred_fclc.opt,"rrBLUP/predict_fclc.opt_trait_1.csv",sep=",",row.names=F,col.names=F)
ability_regression("fclc.opt",pheno_lc,"fclc.opt2lc")

##########################################################################
# GS of BLUP
fw_blup <- read.table("FW-BLUP.csv",head=T,sep=",")[,2]
fw_blup <- matrix(fw_blup,ncol=1)
fcfw_blup <- read.table("FCFW-BLUP.csv",head=T,sep=",")[,2]
fcfw_blup <- matrix(fcfw_blup,ncol=1)

## FW BLUP 
rrBLUP(geno,"blup",fw_blup)
ability_regression("blup",fw_blup,"blup2blup")

## FWFC BLUP
rrBLUP(geno,"fcblup",fcfw_blup)
ability_regression("fcblup",fcfw_blup,"fcblup2fcblup")

##########################################################################

