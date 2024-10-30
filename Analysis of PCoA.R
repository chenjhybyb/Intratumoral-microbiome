
library(ape)
library(vegan)
library(ade4) 
library(tidyverse)



PCoA=function(species,group,cancer,file){
  data <- read.csv(species,header=T,row.names = 1) #log(cpm+1)
  data <- t(data)
  data <- cbind(RNA.id = rownames(data), data)
  rownames(data) <- NULL
  
  group <- read.csv(group,header = T)
 
  sample <- group[,c("RNA.id","group")]
  data2 <- merge(data,sample,by.x="RNA.id",by.y="RNA.id")
  

  data3 <- data2 %>% 
    unite("group_sample", c(group,RNA.id)) %>%
    column_to_rownames("group_sample")
  
df.dist = vegdist(data3,method='euclidean')
pcoa =  dudi.pco(df.dist,
                 scannf = F,   
                 nf=2)        
data = pcoa$li
data <- rownames_to_column(data)
data <- separate(data, rowname, into = c("group", "id.RNA"), sep = "_", remove = TRUE)
adonis_result <- adonis2(df.dist ~ data$group, method = "euclidean")
p_value <- adonis_result$`Pr(>F)`[1]
# 绘图
p <- ggplot(data,aes(x = A1,
                      y = A2,
                      
))+
  geom_point()+
  theme_classic()+
  
  geom_point(aes(color = group),shape = 19,size = 3)+
  scale_color_manual(values = c("#74BBFD","#F24957"),
                     limits = c("R","NR"),name="Group")+
  scale_fill_manual(values = c("#74BBFD","#F24957"),
                    limits = c("R","NR"))+
  
  stat_ellipse(aes(x=A1,    # 添加置信区间圈
                   y=A2,color=group
  ),
  geom = "polygon",
  level = 0.95,
  alpha=0)+
  
  labs(  
    x = paste0("PCoA1 (",as.character(round(pcoa$eig[1] / sum(pcoa$eig) * 100,2)),"%)"),
    y = paste0("PCoA2 (",as.character(round(pcoa$eig[2] / sum(pcoa$eig) * 100,2)),"%)"),title = cancer
  ) + annotate("text", x = Inf, y = Inf, 
               label = paste0("P = ", p_value), size = 6, hjust = 1, vjust = 1)+
  theme(legend.position = "none",axis.title.x = element_text(size = 14),
        axis.title.y =element_text(size = 14),title = element_text(size=14) )

ggsave(
  filename = file, 
  width = 19,     
  height = 4,            
  units = "in",          
  dpi = 300              
)
}

