### BLUP estimation for phenotypic traits across multiple environments
### Variance components of G, E, and GxE
### Broad-sense heritability estimation across multiple environments
### Contact: tonghao0605@gmail.com

#dir <- ("/winmounts/tong/winhome/mpidir/8.AlleleFlux/00_GxE/")
#dir <- ("H:/mpidir/8.AlleleFlux/00_GxE/")
dir <- ("D:/05_fluxGxE/")
setwd(dir)
rm(list=ls())

library(lme4)

#######################################################################
# data set 
dataopt <- read.table("opt/AraCore_allacc_opt.csv",sep=",",header=F)
dataln <- read.table("LN/AraCore_allacc_LN.csv",sep=",",header=F)
datalc <- read.table("LC/AraCore_allacc_LC.csv",sep=",",header=F)
databiom <- read.table("LN/araid67_final_new_env.csv",sep=",",header=F)

#######################################################################
### input data format: row for each line, column for each trait

dataall <- rbind(t(dataopt),t(dataln),t(datalc))
dataall <- cbind(dataall,c(databiom[,2],databiom[,3],databiom[,4]))

n <- ncol(dataopt) #lines number
p <- ncol(dataall) #phenotype number
m <- 3 #environment number

### prepare input data format 
lineid <- rep(1:n,m)
locid <- c(sort(rep(1:m,n)))

xa <- cbind(lineid,locid,dataall)
colnames(xa) <- c("LINE","LOC",1:337)
write.table(xa,"MLM/flux-all-env.csv",sep=',',quote=F,row.names=F)

##########################################################
### BLUP with G-by-E
n <- 67
line.blup <- c(1:n)
heritability <- NULL

x <- read.table("MLM/flux-all-env.csv",sep=",",header=T)
colnames(x)[1:2] <- c("LINE","LOC")

correct <- c(22,23,263)

for(i in 3:ncol(x)){
  
  print(i-2)
  
  # control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
  #                     check.nobs.vs.nlev = "ignore",
  #                     check.nobs.vs.rankZ = "ignore",
  #                     check.nobs.vs.nRE="ignore")
  # varcomp <- lmer(x[,i]~(1|LINE)+(1|LOC)+(1|LINE:LOC),data=x,control=control) 
  # #isSingular(varcomp, tol = 1e-4)
  
  varcomp <- lmer(x[,i]~(1|LINE)+(1|LOC)+(1|LINE:LOC),data=x,  
                  control=lmerControl(check.nobs.vs.nlev = "ignore",
                  check.nobs.vs.rankZ = "ignore",
                  check.nobs.vs.nRE="ignore"))
  
  var.trans <- lme4::VarCorr(varcomp)
  var <- data.frame(Groups=c('LINE','LOC','GxE','Residual'),
                    Variance=c(as.numeric(var.trans$LINE),as.numeric(var.trans$LOC),as.numeric(var.trans$`LINE:LOC`),attr(var.trans,'sc')^2),check.names=F)
  #residual standard deviation is stored as attribute "sc"
  Gvar<-as.numeric(as.character(var$Variance))[var$Groups%in%'LINE']
  Evar<-as.numeric(as.character(var$Variance))[var$Groups%in%'LOC']
  GxEvar<-as.numeric(as.character(var$Variance))[var$Groups%in%'GxE']
  evar<-as.numeric(as.character(var$Variance))[var$Groups%in%'Residual']
  
  if (i==correct[1]+2 | i==correct[2]+2 | i==correct[3]+2){
    Gvar <- 1E-10
  }
  
  h2 <- c(Gvar,Evar,GxEvar,evar,Gvar/(Gvar+Evar+GxEvar+evar),Evar/(Gvar+Evar+GxEvar+evar),
          GxEvar/(Gvar+Evar+GxEvar+evar),evar/(Gvar+Evar+GxEvar+evar),
          Gvar/(Gvar+Evar+GxEvar/m+evar/m))
  heritability <- rbind(heritability,h2)
  
  f <- fixef(varcomp)
  r <- ranef(varcomp)$LINE
  blup <- f+r
  line.blup <- cbind(line.blup,blup)
}
colnames(line.blup) <- c('line',names(x)[-c(1:2)])
heritability <- cbind(c(colnames(x)[-c(1:2)]),heritability)
colnames(heritability) <- c("flux","G variance","E variance","G-by-E varicane","residual","G%","E%","GxE%","e%","h2")

write.table(line.blup,"MLM/flux-BLUP.csv",row.names=F,sep=",")
write.table(heritability,"MLM/flux-h2.csv",row.names=F,col.names=T,sep=",")

##########################################################
