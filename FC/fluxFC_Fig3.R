### check the fold change between environments on fluxes
### Contact: tonghao0605@gmail.com

#dir <- ("/winmounts/tong/winhome/mpidir/8.AlleleFlux/00_GxE/FC")
dir <- ("D:/05_fluxGxE/FC/")
setwd(dir)
rm(list=ls())

library(tidyverse)
library(ggpubr)
library(ggsci)
library(VennDiagram)

#######################################################################
# data set 
dataopt <- read.table("../opt/AraCore_allacc_opt.csv",sep=",",header=F)
dataln <- read.table("../LN/AraCore_allacc_LN.csv",sep=",",header=F)
datalc <- read.table("../LC/AraCore_allacc_LC.csv",sep=",",header=F)

fluxname <- read.table("../fluxname.csv",sep=",",header=F)
fluxsubs <- read.table("../fluxsubsys_336_raw_all.txt",sep="\t",header=F)[1:336,1]

#######################################################################
### fold change

#threshold <- 3
threshold <- 8

# LN 
fcln <- dataln/dataopt
write.table(fcln,"fcln_data.csv",sep=",",row.names=F,col.names=F)

#fcln_ave <- rowMeans(fcln)
fcln_max <- apply(fcln,1,max)
fcln_min <- apply(fcln,1,min)

fcln_min_pos <- NULL
for (i in 1:nrow(fcln)){
  fclni <- fcln[i,]
  fclni_pos <- fclni[which(fclni>0)]
  fclni_min_pos <- min(fclni_pos)
  fcln_min_pos <- c(fcln_min_pos,fclni_min_pos)
}

# define significant of threshold fold change
fcln_pos_diff <- which(fcln_max>=threshold)
fcln_neg_diff <- which(fcln_min_pos<=1/threshold)
fcln_dir_diff <- which(fcln_min<0)
fcln_diff_all <- Reduce(union,list(fcln_pos_diff,fcln_neg_diff,fcln_dir_diff))

fcln_diff_direction <- rep("d",length(fcln_diff_all))
fcln_diff_direction[match(fcln_pos_diff,fcln_diff_all)] <- "+"
fcln_diff_direction[match(fcln_neg_diff,fcln_diff_all)] <- "-"
fcln_diff_direction[match(intersect(fcln_pos_diff,fcln_neg_diff),fcln_diff_all)] <- "+,-"
fcln_diff_direction[match(intersect(fcln_pos_diff,fcln_dir_diff),fcln_diff_all)] <- "+,d"
fcln_diff_direction[match(intersect(fcln_neg_diff,fcln_dir_diff),fcln_diff_all)] <- "-,d"
fcln_diff_direction[match(Reduce(intersect,list(fcln_pos_diff,fcln_neg_diff,fcln_dir_diff)),fcln_diff_all)] <- "+,-,d"

fcln_diff_all_f <- cbind(fcln_diff_all,fluxname[fcln_diff_all,],fluxsubs[fcln_diff_all],
                         fcln_max[fcln_diff_all],fcln_min_pos[fcln_diff_all],fcln_min[fcln_diff_all],
                         fcln_diff_direction)
write.table(fcln_diff_all_f,paste0("fcln_summary_",threshold,"fold.csv"),sep=",",row.names=F,col.names=F)

#######################################################################
# LC
fclc <- datalc/dataopt
write.table(fclc,"fclc_data.csv",sep=",",row.names=F,col.names=F)

#fclc_ave <- rowMeans(fclc)
fclc_max <- apply(fclc,1,max)
fclc_min <- apply(fclc,1,min)

fclc_min_pos <- NULL
for (i in 1:nrow(fclc)){
  fclci <- fclc[i,]
  fclci_pos <- fclci[which(fclci>0)]
  fclci_min_pos <- min(fclci_pos)
  fclc_min_pos <- c(fclc_min_pos,fclci_min_pos)
}

# define significant of threshold fold change
fclc_pos_diff <- which(fclc_max>=threshold)
fclc_neg_diff <- which(fclc_min_pos<=1/threshold)
fclc_dir_diff <- which(fclc_min<0)
fclc_diff_all <- Reduce(union,list(fclc_pos_diff,fclc_neg_diff,fclc_dir_diff))

fclc_diff_direction <- rep("d",length(fclc_diff_all))
fclc_diff_direction[match(fclc_pos_diff,fclc_diff_all)] <- "+"
fclc_diff_direction[match(fclc_neg_diff,fclc_diff_all)] <- "-"
fclc_diff_direction[match(intersect(fclc_pos_diff,fclc_neg_diff),fclc_diff_all)] <- "+,-"
fclc_diff_direction[match(intersect(fclc_pos_diff,fclc_dir_diff),fclc_diff_all)] <- "+,d"
fclc_diff_direction[match(intersect(fclc_neg_diff,fclc_dir_diff),fclc_diff_all)] <- "-,d"
fclc_diff_direction[match(Reduce(intersect,list(fclc_pos_diff,fclc_neg_diff,fclc_dir_diff)),fclc_diff_all)] <- "+,-,d"

fclc_diff_all_f <- cbind(fclc_diff_all,fluxname[fclc_diff_all,],fluxsubs[fclc_diff_all],
                         fclc_max[fclc_diff_all],fclc_min_pos[fclc_diff_all],fclc_min[fclc_diff_all],
                         fclc_diff_direction)
write.table(fclc_diff_all_f,paste0("fclc_summary_",threshold,"fold.csv"),sep=",",row.names=F,col.names=F)

#######################################################################
# overlap between LN and LC
fclnlc_overlap <- intersect(fcln_diff_all,fclc_diff_all)
fclnlc_overlap_f <- cbind(fclnlc_overlap,fluxname[fclnlc_overlap,],fluxsubs[fclnlc_overlap],
                       fcln_diff_all_f[match(fclnlc_overlap,fcln_diff_all),4:7],
                       fclc_diff_all_f[match(fclnlc_overlap,fclc_diff_all),4:7])
write.table(fclnlc_overlap_f,paste0("overlaps_b2_",threshold,"fold.csv"),sep=",",row.names=F,col.names=F)

fcln_only <- setdiff(fcln_diff_all,fclnlc_overlap)
write.table(fcln_only,paste0("fcln_only_b2_",threshold,"fold.csv"),sep=",",row.names=F,col.names=F)

fclc_only <- setdiff(fclc_diff_all,fclnlc_overlap)
write.table(fclc_only,paste0("fclc_only_b2_",threshold,"fold.csv"),sep=",",row.names=F,col.names=F)

#######################################################################
# overlap between LN+- and LC+-
fcln_negdir_diff <- union(fcln_dir_diff,fcln_neg_diff)
fclc_negdir_diff <- union(fclc_dir_diff,fclc_neg_diff)
fclnlc_overlap_pn <- Reduce(intersect,list(fcln_pos_diff,fcln_negdir_diff,fclc_pos_diff,fclc_negdir_diff))

fclnlc_overlap_pn_f <- cbind(fclnlc_overlap_pn,fluxname[fclnlc_overlap_pn,],fluxsubs[fclnlc_overlap_pn],
                          fcln_diff_all_f[match(fclnlc_overlap_pn,fcln_diff_all),4:7],
                          fclc_diff_all_f[match(fclnlc_overlap_pn,fclc_diff_all),4:7])
write.table(fclnlc_overlap_pn_f,paste0("overlaps_b4_",threshold,"fold.csv"),sep=",",row.names=F,col.names=F)
#######################################################################

#######################################################################
# Figure 3a
#######################################################################

##################################
# data set 
dataopt <- read.table("../opt/AraCore_allacc_opt.csv",sep=",",header=F)
dataln <- read.table("../LN/AraCore_allacc_LN.csv",sep=",",header=F)
datalc <- read.table("../LC/AraCore_allacc_LC.csv",sep=",",header=F)

fcln <- dataln/dataopt
fclc <- datalc/dataopt
datafc <- cbind(fcln,fclc)

sub <- read.table("../fluxsubsys_336.csv",sep=",",header=F)
subb <- as.character(sub[,2])
subu <- as.character(unique(sub[,2]))
subuall <- sort(subu)
subfinal <- subuall[2:62]

##################################
# for each subsystem
dataall <- NULL
subsig <- c(6,8,17)
#subsig <- c(1:61)
for (i in subsig){
  subidi <- which(subb==subfinal[i])
  subidii <- sub[subidi,1]
  nsubid <- length(subidi)
  
  datai <- colMeans(datafc[unique(subidii),])
  dataii <- cbind(subfinal[i],datai,c(rep("LN",67),rep("LC",67)))
  dataall <- rbind(dataall,dataii)

}

dataall <- as.data.frame(dataall)
dataall[,2] <- as.numeric(dataall[,2])
colnames(dataall) <- c("subsys","fc","cond")

##################################
## plot 
pheno_fcdata <- dataall[1:134,]
pheno_fcdata <- dataall[135:268,]
pheno_fcdata <- dataall[269:402,]

myCol <- pal_npg("nrc",alpha=1)(10)
phist <- gghistogram(pheno_fcdata, x = "fc", bins = 20,
                     color = "cond", fill = "cond",
                     palette = myCol[c(5,7)], 
                     xlab ="Fold change", ylab = "Number of accession"
) + 
  # scale_y_continuous(expand = expansion(mult = c(0, 0)), position = "left",
  #                    breaks = get_breaks(by = 5, from = 0)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), position = "right",
                     breaks = get_breaks(by = 5, from = 0)) + 
  theme(plot.margin=unit(c(5,8,5,5),"mm"))
phist <- ggpar(phist,font.xtickslab = c(30),font.x = c(30),
               font.ytickslab = c(30),font.y = c(30),
               legend="none", legend.title="",font.legend=c(17))

pdensity <- ggdensity(pheno_fcdata, x = "fc", y ="..count..",
                      color = "cond", size =1.5,
                      palette = myCol[c(5,7)], alpha = 0,
                      xlab ="", ylab = "Density"
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis") + rremove("xlab") + rremove("x.text") + rremove("x.ticks") +
  rremove("y.axis") + rremove("ylab") + rremove("y.text") + rremove("y.ticks") +
  rremove("legend")

aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
pp <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
pp

ggsave("Fig3a-asp.tiff",pp,dpi=300,width=800*300/72,height=400*300/72,units="px",
       device = grDevices::tiff)

data_save <- dataall
colnames(data_save) <- c("Subsys","fold change","condition")
write.table(data_save,"Fig3a.csv",sep=",",row.names=F,col.names=T)

#######################################################################
#######################################################################
#######################################################################
#######################################################################
# Figure 3a-0
#######################################################################

fcln_max_fig <- fcln_max
fcln_min_fig <- fcln_min

#fcln_max_fig[which(fcln_max<0)] <- 2048
#fcln_min_fig[which(fcln_min<0)] <- 1/2048

# fcln_max_log <- log2(fcln_max_fig)
# fcln_min_log <- log2(fcln_min_fig)

corr_all <- NULL
for (i in 1:336){
  cori <- cor(as.numeric(fcln[i,]),as.numeric(fcln[336,]),method = "pearson")
  corr_all <- c(corr_all,cori)
}

group_max <- rep("max",336)
group_min <- rep("min",336)
group_max[which(fcln_max<0)] <- "reverse"
group_min[which(fcln_min<0)] <- "reverse"
dataall <- data.frame(fc=c(fcln_max_fig,fcln_min_fig),r2=rep(corr_all,2),group=c(group_max,group_min))
p <- ggscatter(dataall, x="fc", y="r2",
          color = "group",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          xlab ="Fold change", ylab = "Correlation to fresh weight")+
geom_vline(xintercept=c(-3,3),linetype = 2)
pp <- ggpar(p, legend="top", legend.title="")

ggsave("Fig3a-0.tiff",pp,dpi=300,width=500*300/72,height=300*300/72,units="px")

#######################################################################
# Figure 3a-00
#######################################################################

fclc_max_fig <- fclc_max
fclc_min_fig <- fclc_min

# fclc_max_fig[which(fclc_max<0)] <- 2048
# fclc_min_fig[which(fclc_min<0)] <- 1/2048

# fclc_max_log <- log2(fclc_max_fig)
# fclc_min_log <- log2(fclc_min_fig)

corr_all <- NULL
for (i in 1:336){
  cori <- cor(as.numeric(fclc[i,]),as.numeric(fclc[336,]),method = "pearson")
  corr_all <- c(corr_all,cori)
}

# pv_all <- NULL
# for (i in 1:336){
#   pvi <- t.test(as.numeric(dataopt[i,]),as.numeric(datalc[i,]),paired = T)
#   pv_all <- c(pv_all,-log10(pvi$p.value))
# }

group_max <- rep("max",336)
group_min <- rep("min",336)
group_max[which(fclc_max<0)] <- "reverse"
group_min[which(fclc_min<0)] <- "reverse"
dataall <- data.frame(fc=c(fclc_max_fig,fclc_min_fig),pv=rep(corr_all,2),group=c(group_max,group_min))
p <- ggscatter(dataall, x="fc", y="pv",
               color = "group",
               palette = c("#00AFBB", "#E7B800", "#FC4E07"),
               xlab ="Fold change", ylab = "Correlation to fresh weight")+
geom_vline(xintercept=c(-3,3),linetype = 2)
pp <- ggpar(p, legend="top", legend.title="")

ggsave("Fig3a-00.tiff",pp,dpi=300,width=500*300/72,height=300*300/72,units="px")

#######################################################################
#######################################################################
#######################################################################
#######################################################################


#######################################################################
# Figure 3b
#######################################################################

#######################################################################
# data set 
fclndata <- rep(0,336)
fclndata[fcln_diff_all] <- 1
fclcdata <- rep(0,336)
fclcdata[fclc_diff_all] <- 1

sub <- read.table("../fluxsubsys_336.csv",sep=",",header=F)
subb <- as.character(sub[,2])
subu <- as.character(unique(sub[,2]))
subuall <- sort(subu)
subfinal <- subuall[2:62]

#######################################################################
# for each subsystem
#subidall <- NULL
#nsubidall <- NULL
dataall <- matrix(0,61,2)
for (i in 1:61){
  subidi <- which(subb==subfinal[i])
  subidii <- sub[subidi,1]
  nsubid <- length(subidi)
  
  dataall[i,1] <- length(which(fclndata[unique(subidii)]==1))/length(unique(subidii))
  dataall[i,2] <- length(which(fclcdata[unique(subidii)]==1))/length(unique(subidii))
  
  #subidall <- c(subidall,subidii)
  #nsubidall <- c(nsubidall,nsubid)
}

colnames(dataall) <- c("LN","LC")
dataall <- cbind(name=subfinal,dataall)
dataall <- as.data.frame(dataall)
dataall <- dataall[union(which(dataall[,2]!=0),which(dataall[,3]!=0)),]
data_f <- gather(dataall,"FC.rate","value",c(2:3))
data_f[,3] <- round(as.numeric(data_f[,3]),2)
data_f[which(data_f[,3]==0),3] <- 0.008

#######################################################################
# stacked bar plot
#library(RColorBrewer)
#library(colorspace)
#myCol <- brewer.pal(4, "Pastel2")
#myCol <- darken(myCol, 0.15)
myCol <- pal_npg("nrc",alpha=1)(10)
#myCol <- darken(myCol, 0.05)

p <- ggbarplot(data_f, x = "name", y = "value",
          color = "white", 
          fill = "FC.rate",
          #palette = "npg", 
          palette = myCol[c(5,7)],
          #sort.val = "desc", 
          #sort.by.groups = FALSE,
          position = position_dodge(0.8),
          xlab ="", ylab = "Proportion of reactions with\n FC > 8 or FC < 1/8")+
#rotate_x_text(40)+
theme(plot.margin=unit(c(2,2,2,12),"mm"))
pp <- ggpar(p, font.xtickslab = c(12),
            xtickslab.rt = 40,
            legend="right", legend.title="")
pp
#### save eps in export!!! 
#ggsave("Fig1a.eps",pp,width=1230,height=380,limitsize = FALSE)
#ggsave("Fig1a72.tiff",pp,dpi=72,width=1230,height=380,units="px")
ggsave("Fig3b.tiff",pp,dpi=300,width=600*300/72,height=300*300/72,units="px")

#######################################################################
# Figure 3c
#######################################################################

dataall <- list(ln=as.character(fcln_diff_all),lc=as.character(fclc_diff_all))

# Prepare palette of colors 
# myCol <- pal_npg("nrc",alpha=0.3)(10)
# myCol <- myCol[c(1,3)]
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

# Chart
venn.diagram(
  x = dataall,
  category.names = c("LN","LC"),
  filename = 'Fig3c.tiff',
  output=TRUE,
  
  # Output features
  imagetype="tiff",
  height = 480, 
  width = 480, 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = 0.5,
  #fontface = "bold",
  fontfamily = "sans",
  margin = 0.1,
  
  # Set names
  cat.cex = 0.5,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans"
  #rotation = 1
)

p <- venn.diagram(
  x = dataall,
  category.names = c("LN","LC"),
  filename = NULL,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = 3,
  #fontface = "bold",
  fontfamily = "Arial",
  margin = 0.1,
  
  # Set names
  cat.cex = 3,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "Arial"
  #rotation = 1
)
#grid.newpage()
#grid::grid.draw(p)

ggsave("Fig3c.eps", p, device=cairo_ps)

#######################################################################
# Figure 3d
#######################################################################

dataall <- list(lnpos=as.character(fcln_pos_diff),lnneg=as.character(fcln_negdir_diff),
                lcpos=as.character(fclc_pos_diff),lcneg=as.character(fclc_negdir_diff))

# Prepare palette of colors 
# myCol <- pal_npg("nrc",alpha=0.3)(10)
# myCol <- myCol[c(1,3)]
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

# Chart
venn.diagram(
  x = dataall,
  category.names = c("LN+","LN-","LC+","LC-"),
  filename = 'Fig3d.tiff',
  output=TRUE,
  
  # Output features
  imagetype="tiff",
  height = 480, 
  width = 480, 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:4],
  
  # Numbers
  cex = 0.5,
  #fontface = "bold",
  fontfamily = "sans",
  margin = 0.1,
  
  # Set names
  cat.cex = 0.5,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans"
  #rotation = 1
)

p <- venn.diagram(
  x = dataall,
  category.names = c("LN+","LN-","LC+","LC-"),
  filename = NULL,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:4],
  
  # Numbers
  cex = 3,
  #fontface = "bold",
  fontfamily = "Arial",
  margin = 0.1,
  
  # Set names
  cat.cex = 3,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "Arial"
  #rotation = 1
)
#grid.newpage()
#grid::grid.draw(p)

ggsave("Fig3d.eps", p, device=cairo_ps)

#######################################################################

#######################################################################
# Figure 2-less
#######################################################################

edata <- read.table("enrichment_b2_10fold.csv",sep=",",header=T)
edata$Adjusted.P.value <- -log10(edata$Adjusted.P.value)
edata <- data.frame(edata[-7,])

# library(RColorBrewer)
# myCol <- brewer.pal(4, "Pastel2")
# myCol <- darken(myCol, 0.15)
myCol <- pal_npg("nrc",alpha=1)(10)

p <- ggdotchart(edata, x = "Group", y = "Adjusted.P.value",
                color = "Enriched.set.size",                                # Color by groups
                #palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
                sorting = "descending",                       # Sort value in descending order
                add = "segments",                             # Add segments from y = 0 to dots
                rotate = TRUE,                                # Rotate vertically
                #group = "cyl",                                # Order by groups
                dot.size = 6,                                 # Large dot size
                #label = round(dfm$mpg),                        # Add mpg values as dot labels
                font.label = list(color = "white", size = 9,  # Adjust label parameter
                                  vjust = 0.5),
                xlab ="", ylab = "-log(FDR)") +        
  gradient_color(c(myCol[7],myCol[5]))

pp <- ggpar(p, font.xtickslab = c(13),font.y = c(17),
            legend="top", legend.title="Enriched set size",font.legend=c(9))
pp
#### save eps in export!!! 
#ggsave("Fig1a.eps",pp,width=1230,height=380,limitsize = FALSE)
#ggsave("Fig1a72.tiff",pp,dpi=72,width=1230,height=380,units="px")
ggsave("Fig2-less.tiff",pp,dpi=300,width=600*300/72,height=380*300/72,units="px")

