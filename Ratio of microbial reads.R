
library(tidyverse)
library(plyr)

percentage=function(data1,data2,outfile){
  
  microbes <- read.csv(data1,header = T,row.names = 1)
 
  counts_sum <- data.frame(apply(microbes,MARGIN =2,sum,na.rm = TRUE))
  colnames(counts_sum) <- 'Sum'
  counts_sum$RNA.id <- row.names(counts_sum)
  rownames(counts_sum) <- NULL
  
  host <- read.table(data2)
  colnames(host)[1] <- "RNA.id"
  data <- merge(host,counts_sum,by.x = "RNA.id",by.y = "RNA.id")
  data$percentage=data$Sum/data$V2
  data <- data[,c(1,4)]
  write.csv(data,outfile,row.names = F)
  
}
percentage("PRJNA557841_phylum.csv","PRJNA557841_count.txt","PRJNA557841_percentage.csv")
percentage("EGAD00001004183_phylum.csv","EGAD00001004183_read_count.txt","EGAD00001004183_percentage.csv")
percentage("PRJEB23709_phylum.csv","PRJEB23709_count.txt","PRJEB23709_percentage.csv")
percentage("PRJEB25780_phylum.csv","PRJEB25780_count.txt","PRJEB25780_percentage.csv")
percentage("PRJNA312948_phylum.csv","PRJNA312948_count.txt","PRJNA312948_percentage.csv")
percentage("PRJNA356761_phylum.csv","PRJNA356761_count.txt","PRJNA356761_percentage.csv")
percentage("PRJNA694637_phylum.csv","PRJNA694637_count.txt","PRJNA694637_percentage.csv")

NSCLC <- read.csv("PRJNA557841_percentage.csv",header = T)
NSCLC$type="NSCLC"

RCC <- read.csv("EGAD00001004183__percentage.csv",header = T)
RCC$type = "RCC"
RCC_group <- read.csv("group/EGAD00001004183_allPD-1.csv",header = T)
sample <- RCC_group$RNA.id
RCC <- RCC[RCC$RNA.id %in% sample,]

Gide <- read.csv("PRJEB23709_percentage.csv",header = T)
Gide$type = "Melanoma_Gide"
Gide_group <- read.csv("group/Gide.csv",header = T)
sample <- Gide_group$RNA.id
Gide <- Gide[Gide$RNA.id %in% sample,]

melanoma_group <- read.csv("melanoma_merge/pre/pre_metadata.csv",header = T)
sample <- melanoma_group$id_RNA
Riaz <- read.csv("PRJNA356761_percentage.csv",header = T)
Riaz$type = "Melanoma_Riaz"
Riaz <- Riaz[Riaz$RNA.id %in% sample,]

Hugo <- read.csv("PRJNA312948_percentage.csv",header = T)
Hugo$type = "Melanoma_Hugo"
Hugo <- Hugo[Hugo$RNA.id %in% sample,]


gc <- read.csv("PRJEB25780_percentage.csv",header = T)
gc$type = "Gastric cancer"
gc_group <- read.csv("group/PRJEB25780.csv",header = T)
sample <- gc_group$RNA.id
gc <- gc[gc$RNA.id %in% sample,]

EC <- read.csv("PRJNA694637.csv",header = T)
EC$type = "Esophageal adenocarcinoma"
EC_group <- read.csv("group/PRJNA694637.csv",header = T)
sample <- EC_group$id_RNA
EC <- EC[EC$RNA.id %in% sample,]

data <- rbind(gc,EC,NSCLC,Gide,Hugo,Riaz,RCC)
data$type <- factor(data$type,levels=c("rEAC","mGC","NSCLC",
                                       "Melanoma_Gide","Melanoma_Riaz","Melanoma_Hugo","RCC"))

p <- ggplot(data,aes(x = type, y = percentage_log, fill = type)) +
  geom_dotplot(binaxis = "y",   
               binwidth = 0.12,    
               stackdir = "center",stroke = NA) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", fatten = 2, width = 0.5, color = "black", alpha = 0.50) +
  theme_classic2() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),axis.text.y = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=16,face = "bold")) +
  scale_fill_manual(values = c("#FEBA07","#BB97C5","#EA8958","#2983BB","#97C24E","#207F4C","#EE4863")) +
  xlab("") +
  ylab('Fractional bacteria reads (log10)')
p
