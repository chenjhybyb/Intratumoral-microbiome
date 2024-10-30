
library(Hmisc)

###Calculate immune score
immune_score=function(datain,gene_list,signature_score,out){
  data <- read.csv(datain,header=T)# RNA-seq data

  colnames(data)[1]<-"gene"
  immune_gene<-read.csv(gene_list,header = T)
  colnames(data)[1]<-"gene"
  data2<-merge(immune_gene,data,by.x="gene",by.y="gene")
  data2 <- data2[apply(data2, 1, function(row) !all(row == 0)), ]
  rownames(data2) <- data2[,1]
  data2 <- data2[,-1]
  
  exp <- t(data2)
  exp_data<-data.frame(exp)
  exp_data$immune_score<-apply(exp_data[,seq(1,ncol(exp_data))],1,mean,na.rm=T)
  exp_pre_score<-exp_data[,ncol(exp_data),drop=FALSE]
  colnames(exp_pre_score)[1] <- signature_score

  write.csv(exp_pre_score,out)
}

##picture
library(ggpubr)
library(jjPlot)
library(reshape2)
library(tidyverse)
group <- read.csv(group,header = T)
colnames(MHC)[1] <- "RNA.id"
colnames(IFN)[1] <- "RNA.id"
colnames(IC)[1] <- "RNA.id"
colnames(CYT)[1] <- "RNA.id"
data <- merge(IC, IFN,by.x="RNA.id",by.y="RNA.id")
data <- merge(data, CYT,by.x="RNA.id",by.y="RNA.id")
data <- merge(data, MHC,by.x="RNA.id",by.y="RNA.id")
data <- merge(data,group,by.x="RNA.id",by.y="RNA.id")

data <- pivot_longer(data = data,
                     cols = 2:7,
                     names_to = "Immune_signatures",
                     values_to = "Score")

data$Immune_signatures <- factor(data$Immune_signatures,levels = c("PDL1","PD1","CTLA4",
                                                                   "CYT","IFNγ",
                                                                   "MHCI"))

data$type <- ifelse(data$group == "NR", "left", "right")

p <-  ggplot(data,aes(x = Immune_signatures  , y = Score)) +
  geom_jjviolin(aes(fill = group ,type= type,color=group),
                trim = F,
                width = 0.2,
                shift = -0.025,
                position = position_dodge(width = 0.72),alpha=0.8) +
  scale_fill_manual(values = c('#F24957','#74BBFD'),name="Group") +
  scale_color_manual(values = c('#F24957','#74BBFD'),name="Group")+
  theme_classic()+theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),axis.text.x = element_text(size = 15,face = "bold"),
                        axis.title.x = element_blank() ,legend.position = "none",axis.title.y = element_text(size=14,face = "bold"),
                        axis.text.y = element_text(size=12),
                        title = element_text(size=15) )       +                          
  
  stat_compare_means(aes(group=group), label="p.format",method="wilcox.test",label.y = max(data$Score) + 1.4)
p
