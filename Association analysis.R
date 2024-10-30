

###Logistic regression analysis

##Response
response<-function(micro,clindata,out){
  
  # micro<-read.csv(micro,header=T)
  # clin<-read.csv("1_clinical/Mariathasan2018_clinical.csv",header=T)
  
  micro<-read.csv(micro,header=T)
  clin<-read.csv(clindata,header=T)
  colnames(micro)[1]<-"id_RNA"
  clin<-subset(clin,is.na(id_RNA)==F,drop=T)
  
  microbe_ID<-colnames(micro)[-1]
  nmicrobe<-ncol(micro)-1
  
  pvalue<-c()
  or<-c()
  orl<-c()   # 95% CI left 
  orr<-c()   # 95% CI right
  se<-c()
  
  #i=4
  
  for(i in 1:nmicrobe){
    microbe<-micro[,c(1,(i+1))]
    colnames(microbe)[2]<-"PC"
    cor<-merge(clin,microbe,by.clin="id_RNA",by.microbe="id_RNA")
    
    # ???????ͼ???????��
    if(nrow(cor)>=15 & length(table(cor$response))>1){
      
      
      value<-cor$PC
   
      
      qt<-quantile(value, probs = c(0.25,0.5,0.75),na.rm=T)
      n<-nrow(cor)
      
      for(j in 1:n){
        if(is.na(value[j])==T){ cor$Qt[j]=NA}
        else if(value[j]<=qt[1]){ cor$Qt[j]=1}
        else if (value[j]>qt[1] & value[j]<=qt[2]){cor$Qt[j]=2}
        else if (value[j]>qt[2] & value[j]<=qt[3]){cor$Qt[j]=3}
        else if (value[j]>qt[3]){cor$Qt[j]=4}
      }
      
      a<-glm(response~Qt, family=binomial(logit), data=cor)
      ci<-confint(a)
      ci<-unlist(ci)
      
      b<-unlist(summary(a))
      pvalue[i]<-b$coefficients8
      or[i]<-exp(b$coefficients2) 
      
      orl[i]<-exp(ci[2])
      orr[i]<-exp(ci[4])
      se[i]<-b$coefficients4
    }
    else {
      pvalue[i]<-NA
      or[i]<-NA
      
      orl[i]<-NA
      orr[i]<-NA
      se[i]<-NA
    }
  }
  pvalue<-as.matrix(pvalue,nmicrobe,1)
  or<-as.matrix(or,nmicrobe,1)
  orl<-as.matrix(orl,microbe,1)
  orr<-as.matrix(orr,nmicrobe,1)
  se<-as.matrix(se,nmicrobe,1)
  
  re<-cbind(pvalue,or,orl,orr,se)
  
  rownames(re)<-microbe_ID
  colnames(re)<-c("pvalue","or","or_left","or_right","se")
  re<-data.frame(re)
  write.csv(re,out)
}

##PFS
pfs<-function(micro,clindata,out){
  
  
  micro<-read.csv(micro,header=T)
  clin<-read.csv(clindata,header=T)
  colnames(micro)[1]<-"id_RNA"
  clin<-subset(clin,is.na(id_RNA)==F,drop=T)
  
  microbe_ID<-colnames(micro)[-1]
  nmicrobe<-ncol(micro)-1
  
  pvalue<-c()
  or<-c()
  orl<-c()   # 95% CI left 
  orr<-c()   # 95% CI right
  se<-c()
  
  for(i in 1:nmicrobe){
    microbe<-micro[,c(1,(i+1))]
    colnames(microbe)[2]<-"PC"
    cor<-merge(clin,microbe,by.clin="id_RNA",by.microbe="id_RNA")
    

    if(nrow(cor)>=15){
      
      
      value<-cor$PC
     
      
      qt<-quantile(value, probs = c(0.25,0.5,0.75),na.rm=T)
      n<-nrow(cor)
      
      for(j in 1:n){
        if(is.na(value[j])==T){ cor$Qt[j]=NA}
        else if(value[j]<=qt[1]){ cor$Qt[j]=1}
        else if (value[j]>qt[1] & value[j]<=qt[2]){cor$Qt[j]=2}
        else if (value[j]>qt[2] & value[j]<=qt[3]){cor$Qt[j]=3}
        else if (value[j]>qt[3]){cor$Qt[j]=4}
      }
      
      a<-summary(coxph(Surv(pfs, pfs_event)~Qt, cor,na.action=na.exclude))
      a<-as.vector(unlist(a))
      pvalue[i]<-a$coefficients5
      HR[i]<-a$coefficients2
      HR_left[i]<-a$conf.int3
      HR_right[i]<-a$conf.int4
      se[i]<-a$coefficients3
      
    }
    else {
      pvalue[i]<-NA
      HR[i]<-NA
      HR_left[i]<-NA
      HR_right[i]<-NA
      se[i]<-NA
    }
  }
  pvalue<-as.matrix(pvalue,nmicrobe,1)
  HR<-as.matrix(HR,nmicrobe,1)
  HR_left<-as.matrix(HR_left,nmicrobe,1)
  HR_right<-as.matrix(HR_right,nmicrobe,1)
  se<-as.matrix(se,nmicrobe,1)
  re<-data.frame(cbind(pvalue,HR,HR_left,HR_right,se))
  
  rownames(re)<-microbe_ID
  colnames(re)<-c("pvalue","HR","HR_left","HR_right","se")
  re<-data.frame(re)
  
  ######## raw p-value
  write.csv(re,out)
}

##OS
os<-function(micro,clindata,out){
  micro<-read.csv(micro,header=T)
  clin<-read.csv(clindata,header=T)
  colnames(micro)[1]<-"id_RNA"
  clin<-subset(clin,is.na(id_RNA)==F,drop=T)
  
  microbe_ID<-colnames(micro)[-1]
  nmicrobe<-ncol(micro)-1
  
  pvalue<-c()
  or<-c()
  orl<-c()   # 95% CI left 
  orr<-c()   # 95% CI right
  se<-c()
  
  for(i in 1:nmicrobe){
    microbe<-micro[,c(1,(i+1))]
    colnames(microbe)[2]<-"PC"
    cor<-merge(clin,microbe,by.clin="id_RNA",by.microbe="id_RNA")
    
    if(nrow(cor)>=15){
      
   
      value<-cor$PC

      
      qt<-quantile(value, probs = c(0.25,0.5,0.75),na.rm=T)
      n<-nrow(cor)
      
      for(j in 1:n){
        if(is.na(value[j])==T){ cor$Qt[j]=NA}
        else if(value[j]<=qt[1]){ cor$Qt[j]=1}
        else if (value[j]>qt[1] & value[j]<=qt[2]){cor$Qt[j]=2}
        else if (value[j]>qt[2] & value[j]<=qt[3]){cor$Qt[j]=3}
        else if (value[j]>qt[3]){cor$Qt[j]=4}
      }
      
      a<-summary(coxph(Surv(os, os_event)~Qt, cor,na.action=na.exclude))
      a<-as.vector(unlist(a))
      pvalue[i]<-a$coefficients5
      HR[i]<-a$coefficients2
      HR_left[i]<-a$conf.int3
      HR_right[i]<-a$conf.int4
      se[i]<-a$coefficients3
      
    }
    else {
      pvalue[i]<-NA
      HR[i]<-NA
      HR_left[i]<-NA
      HR_right[i]<-NA
      se[i]<-NA
    }
  }
  pvalue<-as.matrix(pvalue,nmicrobe,1)
  HR<-as.matrix(HR,nmicrobe,1)
  HR_left<-as.matrix(HR_left,nmicrobe,1)
  HR_right<-as.matrix(HR_right,nmicrobe,1)
  se<-as.matrix(se,nmicrobe,1)
  re<-data.frame(cbind(pvalue,HR,HR_left,HR_right,se))
  
  rownames(re)<-microbe_ID
  colnames(re)<-c("pvalue","HR","HR_left","HR_right","se")
  re<-data.frame(re)
  
  ######## raw p-value
  write.csv(re,out)
}

##Spearman correlation analysis
cor.mtest <- function(microbe, clin,outfile){
mat1 <- read.csv(microbe,header = T,row.names = 1)
mat2 <- read.csv(clin,header = T,row.names = 1)
mat2 <- mat2[rownames(mat1),]
mat1 <- mat1[sort(rownames(mat1)),]
mat2 <- mat2[sort(rownames(mat2)),]
mat2 <- mat2[,c(4,6)]
pair_wise_res <- as.data.table(expand.grid(colnames(mat2), colnames(mat1)))
colnames(pair_wise_res) <- c('feature', 'taxa')
pair_wise_res[, c("rho", "rho_p", "tau", "tau_p"):=list(numeric(), numeric(), numeric(), numeric())]
setkey(pair_wise_res, feature, taxa)
for (mat1_idx in 1:ncol(mat1)) {
  for (mat2_idx in 1:ncol(mat2)) {
    cell_name <- colnames(mat2)[mat2_idx]
    microbe_type <- colnames(mat1)[mat1_idx]
    rho_test_res <- cor.test(mat1[,mat1_idx], mat2[,mat2_idx],
                             method="spearman")
    tau_test_res <- cor.test(mat1[,mat1_idx], mat2[,mat2_idx],
                             method='kendall', continuity=T)
    pair_wise_res[.(feature_name, microbe_type), rho_p:=rho_test_res$p.value]
    pair_wise_res[.(feature_name, microbe_type), rho:=rho_test_res$estimate]
    pair_wise_res[.(feature_name, microbe_type), tau_p:=tau_test_res$p.value]
    pair_wise_res[.(feature_name, microbe_type), tau:=tau_test_res$estimate]
  }
}
hist(pair_wise_res$rho_p,breaks = 30)
write.csv(outfile,row.names = F)
}