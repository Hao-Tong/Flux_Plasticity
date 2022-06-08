
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- ("D:/05_fluxGxE/FC_LNLC/rrBLUP")
setwd(dir)

options(digits = 15)

##########################################################################
## load data
model <- c("fcln2fcln","fcln.opt2ln","fclc2fclc","fclc.opt2lc",
           "fcln2fclc","fcln.opt2lc","fclc2fcln","fclc.opt2ln")

outputall <- NULL
for (i in 1:length(model)){
  r2i <- read.table(paste0("summary_",model[i],".csv"),sep=",",header=T)[,2]
  outputall <- cbind(outputall,r2i)
}

outputall <- cbind(c(1:336),outputall)
colnames(outputall) <- c("ID",model)
write.table(outputall,"00_summary_all_sets_group.csv",sep=",",row.names=F,col.names=T)

##########################################################################
