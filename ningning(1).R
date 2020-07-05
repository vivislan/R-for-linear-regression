rm(list = ls())
setwd("/Users/xiaoqi/collaborations/ningning/result/cor_plot")
dat0 = read.delim("/Users/xiaoqi/collaborations/ningning/result/detail.txt")
dat = dat0[,-c(12,14,15,18:20)]
rownames(dat) = dat$姓名


pdf("corr_summary_reg_all.pdf",9,9)
par(mar = c(3.5, 3, 1.6, 1.1), mgp = c(1.9, 0.5, 0),mfrow = c(2,2))

for (lv in 2:7){
  for (ft in colnames(dat)[10:21]){
    # lv = 6
    # ft = colnames(dat)[11]

    pval_mat = get_corr_pval(lv = lv, ft = ft)
    #write.csv(pval_mat,file = paste0("Corr_matrix_",lv,"_",ft),quote = F,eol = "\n",row.names = T,col.names = T)
  }
}

dev.off()


get_corr_pval <- function(lv,ft){
  fn = paste0("/Users/xiaoqi/collaborations/ningning/result/02.OTUanalysis/taxa-bar-plots/taxa-bar-plots.ITS1.sjtu/level-",lv,".csv")
  virus0 = read.csv(fn,row.names = 1)
  
  POS = c(CV = "口",CL = "中")
  for (p in 1:2){
    virus = virus0[virus0$group3==POS[p],]
    row.names(virus) = virus$姓名
    virus = virus[,grep(colnames(virus0),pattern = "k__")]
    virus = virus/rowSums(virus)
    
    
    # ## remove the most redundency
    # virus = virus[,colSums(virus)/sum(virus) < 0.9] 
    # virus = virus/rowSums(virus)
    
    comm = intersect(rownames(dat),rownames(virus))
    for(v in colnames(virus)){
      if (var(virus[comm,v]) > 0){
        
        ## regression remove confounding factors
        regdata = cbind(prop = virus[comm,v],dat[comm,c("年龄","身高","体重","BMI","AFS评分")],ft=dat[comm,ft])
        lm.out = lm(formula = prop ~ ., data = regdata)
        RegPval = summary(lm.out)$coefficients["ft","Pr(>|t|)"]
        
        ## pvalue 和 qvalue
        Ppval = round(cor.test(virus[comm,v],dat[comm,ft],method = "pearson")[["p.value"]],4)
        Spval = round(cor.test(virus[comm,v],dat[comm,ft],method = "spearman")[["p.value"]],4)
        Pearson = round(cor(virus[comm,v],dat[comm,ft],method = "pearson"),4)
        Spearman = round(cor(virus[comm,v],dat[comm,ft],method = "spearman"),4)
        #if(Ppval < 0.05 & Spval < 0.05 ){ #& abs(Pearson) > 0.5 & abs(Spearman) > 0.5
        if(RegPval < 0.05){ #& abs(Pearson) > 0.5 & abs(Spearman) > 0.5
          v.short = tail(strsplit(v,split = "\\.")[[1]],1)
          plot(virus[comm,v],dat[comm,ft],col="#00000050",pch = 19,main = paste0(names(POS)[p]," samples_","level_",lv,"_",ft),xlab = paste0("proportion of '",v,"'"),ylab = paste0("Expression of ",ft),cex.lab = 0.8)
          abline(lm(dat[comm,ft] ~ virus[comm,v]),lwd = 2,lty = 2,col = "gray")
          text(x = (max(virus[comm,v])+min(virus[comm,v]))/2,y = (max(dat[comm,ft])+min(dat[comm,ft]))/2,labels = paste0("Pearson Corr = ",Pearson,"\nPearson Pval = ",Ppval,"\nSpearman Corr = ",Spearman,"\nSpearman Pval = ",Spval,"\n\nRegPval =",RegPval))
          
        }
      }
    }
  }
  # return(pval_mat)
}


######## 0311 #########
fn = paste0("./02.OTUanalysis/taxa-bar-plots/taxa-bar-plots.ITS1.sjtu/level-7.csv")
virus0 = read.csv(fn,row.names = 1)
# virus = virus0[,grep(colnames(virus0),pattern = "k__")]

v = "s__Candida_parapsilosis" #"s__Candida_maltosa"
virus = virus0[virus0$group3=="口",]
row.names(virus) = virus$姓名

commsample = intersect(rownames(virus),rownames(dat))
plot(virus[commsample,grep(pattern = v,x = colnames(virus))],dat[commsample,"IL.6"])



#######################
### pvalues for level 7
fn = paste0("./02.OTUanalysis/taxa-bar-plots/taxa-bar-plots.ITS1.sjtu/level-7.csv")
virus0 = read.csv(fn,row.names = 1)
virus = virus0[,grep(colnames(virus0),pattern = "k__")]
pval_mat_CL = c()
pval_mat_CV = c()

for (ft in colnames(dat)[10:21]){
  pvector = get_corr_mat(lv = 7, ft = ft, p = 1)
  pval_mat_CV = cbind(pval_mat_CV,pvector)
}
rownames(pval_mat_CV) = colnames(virus)
colnames(pval_mat_CV) = colnames(dat)[10:21]

for (ft in colnames(dat)[10:21]){
  pvector = get_corr_mat(lv = 7, ft = ft, p = 2)
  pval_mat_CL = cbind(pval_mat_CL,pvector)
}
rownames(pval_mat_CL) = colnames(virus)
colnames(pval_mat_CL) = colnames(dat)[10:21]

library(gplots)
library(RColorBrewer)
library(pheatmap)
pval_mat_CV=pval_mat_CV[rowSums(is.nan(pval_mat_CV)) == 0,]
pval_mat_CV = pval_mat_CV[rowSums(pval_mat_CV > (-log10(0.05))) >= 1,]

rownames(pval_mat_CV) = 1:16

pval_mat_CV[pval_mat_CV < (-log10(0.05))] = NA
pheatmap(pval_mat_CV,cluster_rows=F,cluster_cols = F,col =  brewer.pal(9, "YlOrRd"))


get_corr_mat <- function(lv,ft,p){
  # fn = paste0("./02.OTUanalysis/taxa-bar-plots/taxa-bar-plots.ITS1.sjtu/level-",lv,".csv")
  # virus0 = read.csv(fn,row.names = 1)
  
  POS = c(CV = "口",CL = "中")
  
  virus = virus0[virus0$group3==POS[p],]
  row.names(virus) = virus$姓名
  virus = virus[,grep(colnames(virus0),pattern = "k__")]
  virus = virus/rowSums(virus)
  
  pvector = c()
  comm = intersect(rownames(dat),rownames(virus))
  for(v in colnames(virus)){
    ## regression remove confounding factors
    regdata = cbind(prop = virus[comm,v],dat[comm,c("年龄","身高","体重","BMI","AFS评分")],ft=dat[comm,ft])
    lm.out = lm(formula = prop ~ ., data = regdata)
    RegPval = summary(lm.out)$coefficients["ft","Pr(>|t|)"]
    RegPval = -log10(RegPval)
    # if(var(virus[comm,v]) > 0){
    #   if (cor(virus[comm,v],dat[comm,ft]) < 0){
    #     RegPval = RegPval*(-1)
    #   }
    # }
    pvector = c(pvector,RegPval)
  }
  return(pvector)
}

