
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- ("D:/05_fluxGxE/FW/rrBLUP")
setwd(dir)

options(digits = 15)

##########################################################################
## load data
model <- c("opt2opt","fcln2fcln","opt2ln","fcln.opt2ln","fcln.rawopt2ln","fclc2fclc","opt2lc","fclc.opt2lc","fclc.rawopt2lc")

outputall <- NULL
for (i in 1:length(model)){
  r2i <- read.table(paste0("summary_",model[i],".csv"),sep=",",header=T)
  outputall <- rbind(outputall,r2i)
}

rownames(outputall) <- model
write.table(outputall,"00_summary_all_sets_group.csv",sep=",",row.names=T,col.names=F)
##########################################################################
