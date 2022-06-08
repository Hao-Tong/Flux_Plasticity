
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
#dir <- ("D:/05_fluxGxE/FC_LNLC/rrBLUP")
#dir <- ("/winmounts/tong/winhome/mpidir/8.AlleleFlux/00_GxE/FC_LNLC/rrBLUP")
dir <- "H:/mpidir/8.AlleleFlux/00_GxE/FC_LNLC/rrBLUP_within"
setwd(dir)

options(digits = 15)

##########################################################################
## LN ##

dataall <- NULL
for (i in 1:336){
  data <- read.table(paste0("predict_fcln.opt_trait_",i,".csv"),sep=",",header=F)
  dataall <- rbind(dataall,data)
}

for (i in 1:50){
  datai <- dataall[,i]
  
  outputall <- NULL
  for (j in 1:336){
    outputall <- cbind(outputall,datai[(67*(j-1)+1):(67*j)])
  }
  
  write.table(outputall,paste0("predict_fcln.opt_rep_",i,".csv"),sep=",",row.names=F,col.names=F)
}

##########################################################################
##########################################################################
## LC ##

dataall <- NULL
for (i in 1:336){
  data <- read.table(paste0("predict_fclc.opt_trait_",i,".csv"),sep=",",header=F)
  dataall <- rbind(dataall,data)
}

for (i in 1:50){
  datai <- dataall[,i]
  
  outputall <- NULL
  for (j in 1:336){
    outputall <- cbind(outputall,datai[(67*(j-1)+1):(67*j)])
  }
  
  write.table(outputall,paste0("predict_fclc.opt_rep_",i,".csv"),sep=",",row.names=F,col.names=F)
}

##########################################################################

##########################################################################
## LN ##

dataall <- NULL
for (i in 1:336){
  data <- read.table(paste0("predict_sign_fcln.opt_trait_",i,".csv"),sep=",",header=F)
  dataall <- rbind(dataall,data)
}

for (i in 1:50){
  datai <- dataall[,i]
  
  outputall <- NULL
  for (j in 1:336){
    outputall <- cbind(outputall,datai[(67*(j-1)+1):(67*j)])
  }
  
  write.table(outputall,paste0("predict_sign_fcln.opt_rep_",i,".csv"),sep=",",row.names=F,col.names=F)
}

##########################################################################
##########################################################################
## LC ##

dataall <- NULL
for (i in 1:336){
  data <- read.table(paste0("predict_sign_fclc.opt_trait_",i,".csv"),sep=",",header=F)
  dataall <- rbind(dataall,data)
}

for (i in 1:50){
  datai <- dataall[,i]
  
  outputall <- NULL
  for (j in 1:336){
    outputall <- cbind(outputall,datai[(67*(j-1)+1):(67*j)])
  }
  
  write.table(outputall,paste0("predict_sign_fclc.opt_rep_",i,".csv"),sep=",",row.names=F,col.names=F)
}

##########################################################################
##########################################################################
# internal cross-validation
##########################################################################
## LN ##

dataall <- NULL
for (i in 1:336){
  data <- read.table(paste0("pearson_LN_within_trait_",i,".csv"),sep=",",header=F)[16,]
  dataall <- rbind(dataall,data)
}

write.table(t(dataall),"pearson_fcln.opt2ln_within.csv",sep=",",row.names=F,col.names=F)

##########################################################################

## LC ##

dataall <- NULL
for (i in 1:336){
  data <- read.table(paste0("pearson_LC_within_trait_",i,".csv"),sep=",",header=F)[16,]
  dataall <- rbind(dataall,data)
}

write.table(t(dataall),"pearson_fclc.opt2lc_within.csv",sep=",",row.names=F,col.names=F)

##########################################################################



