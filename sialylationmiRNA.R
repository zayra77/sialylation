setwd("X:/breast cancer")
.libPaths("D:/R/library")

rm(list=ls())
library(rms)
library(openxlsx)
library(survival)
library(dplyr)
library(survival)
library(survminer)
library(stringr)
library(ggplot2)
library(ggpubr)##add p value
library(rstatix)
library(ggtext)
library(data.table)
library(Rfast)
library(Hmisc)

mRNA<-read.table("X:/tcga data/breast/exp")
miRNA_pre<-read.table("X:/tcga data/breast/mirna")

rawmiRNA<-read.csv("./TCGAmiRNA/TCGA-BRCA_miRNA_matrue_RPM.csv")
rawmiRNA<-rawmiRNA[,-1]
rownames(rawmiRNA)<-rawmiRNA[,1]
rawmiRNA<-rawmiRNA[,-1]
colnames(rawmiRNA)<-substr(colnames(rawmiRNA),0,15)
zerop<-0.65# keep columns with fewer than this threshold of zero values
zeroproportion<-colMeans(rawmiRNA==0)
columns_to_keep<-zeroproportion<zerop

miRNA_less_zero<-rawmiRNA[,columns_to_keep]
miRNA<-miRNA_less_zero+1 #pseudocount
miRNA<-log(miRNA)

library(stringr)
genesymbol<-str_replace_all(rownames(mRNA),"\\|\\d+","")
mRNA<-mRNA[-c(1:30),]
genesymbol<-genesymbol[-c(1:30)]
rownames(mRNA)
#grep("SLC35E2",genesymbol)
sialylationGene<-read.csv("./sialylationgenelist.csv")
sialylationGene<-sialylationGene[,-1]

sialylationrow<-which(genesymbol%in% sialylationGene$Gene.Symbol)

sialylationmRNA<-mRNA[sialylationrow,]
siaSymbol<-str_replace_all(rownames(sialylationmRNA),"\\|\\d+","")
siamRNA<-sialylationmRNA[,intersect(colnames(mRNA),colnames(miRNA))]
siamRNA.df<-data.frame(t(siamRNA))
siamRNA.df<-siamRNA.df+1
siamRNA.df<-log(siamRNA.df)
siamiRNA<-miRNA[,intersect(colnames(mRNA),colnames(miRNA))]
siamiRNA.matrix<-as.matrix(siamiRNA)
siamRNA.matrix<-as.matrix(siamRNA)
# sialylationGene[(!(sialylationGene$Gene.Symbol %in% siaSymbol)),]
cor.mRNA.miRNA<-rcorr(t(siamiRNA.matrix),t(siamRNA.matrix),type="pearson")
cor.mRNA.miRNA$r<-cor.mRNA.miRNA$r[-c(1:nrow(miRNA)),]
cor.mRNA.miRNA$r<-cor.mRNA.miRNA$r[,-c((nrow(miRNA)+1):(nrow(cor.mRNA.miRNA$n)+1))]
cor.mRNA.miRNA$P<-cor.mRNA.miRNA$P[-c(1:nrow(miRNA)),]
cor.mRNA.miRNA$P<-cor.mRNA.miRNA$P[,-c((nrow(miRNA)+1):(nrow(cor.mRNA.miRNA$n)+1))]

cc3<-which(abs(cor.mRNA.miRNA$r)>0.3,arr.ind = TRUE)

cor.mRNA.miRNA$r[cc3]
cor.mRNA.miRNA$P[cc3]
cc3<-as.data.frame(cc3)
rownames(siamiRNA[unique(cc3$col),]) ## 311 sialylation miRNA
cctemp<-which(abs(cor.mRNA.miRNA$r)>0.3 & cor.mRNA.miRNA$P<0.05,arr.ind = TRUE)
cctemp<-as.data.frame(cctemp)
rownames(siamiRNA[unique(cctemp$col),]) ## 910 sialylation miRNA



#####prognosic
survival<-read.table("X:/tcga data/breast/survival",header=TRUE)
miRNA311<-siamiRNA[rownames(siamiRNA[unique(cc3$col),]),]
survival$PatientID<-toupper(survival$PatientID)
survivalunique<-survival %>% distinct()
rownames(survivalunique)<-survivalunique$PatientID
survivalunique<-survivalunique[substr(colnames(miRNA311),0,12),]
miRNA311<-t(miRNA311)
miRNA311<-data.frame(miRNA311)
#miRNA311.log <- miRNA311 %>% mutate_at(c(1:311),log10) 
miRNA311.log<-miRNA311
miRNA311.log<-as.matrix(miRNA311.log)

miRNA311.noinf<-miRNA311.log[,colSums(is.infinite(miRNA311.log))==0]
miRNA311.survival<-cbind(miRNA311.noinf,survivalunique)
mRNA311.survival<-cbind(siamRNA.df,survivalunique)
mRNA311.survival<-within(mRNA311.survival,rm(PatientID))
miRNA311.survival<-within(miRNA311.survival,rm(PatientID))
miRNA311.survival$Survival<-miRNA311.survival$Survival/365
mRNA311.survival$Survival<-mRNA311.survival$Survival/365
p.df<-data.frame(t=numeric(0),OStrain=numeric(0), RFStrain=numeric(0), OStest=numeric(0),RFStest=numeric(0), RFScombine=numeric(0))

for (t in 19:100) {
  
set.seed(20) #77 very good#16 OS ok#50 no 70 no 10 no 20 no 30no 40 good 60 no 80 good 90 no 100 no
Samples<-sample(c(rep(0, 0.8*nrow(miRNA311.survival)), rep(1, 0.2*nrow(miRNA311.survival))))
miRNA311.survival.train<-miRNA311.survival[Samples==0,]
mRNA311.survival.train<-mRNA311.survival[Samples==0,]
miRNA311.survival.test<-miRNA311.survival[Samples==1,]
mRNA311.survival.test<-mRNA311.survival[Samples==1,]
miRNA311.survival.train<-miRNA311.survival.train[colMeans(miRNA311.survival.train==0)<0.98]

pFilter=0.05 #??һ??pֵ??׼????????
outResult=data.frame() #??һ???հ????ݿ򣬺???forѭ????????

sigGenes=c("Death","Survival") #??һ????��??????forѭ???????ã???Ϊ???滹Ҫ?õ?surstat??surtime???????ȷ?????��??
train.column<-ncol(miRNA311.survival.train)-2
for(i in colnames(miRNA311.survival.train[,1:train.column])){ 
  tdcox <- coxph(Surv(Survival, Death) ~ miRNA311.survival.train[,i], data = miRNA311.survival.train)#??ʼ??һѭ??cox????
  tdcoxSummary = summary(tdcox) #summary??????tdcox?ܽᣬ??????????ȡ????
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #??ȡpֵ????????Ҫ?Ǻ?????ȡ????????gene??
  print(pvalue)
  print(i)
  if(pvalue< pFilter){ # ???????ǲ???Ҫ??ȡ???л????????ݣ?ֻ??Ҫ????????gene?Ľ?????????????????pvalue<0.05????ȡ????
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,#?ϲ??У?ʵ?????Ƕ?ѭ???????ĺϲ???ǰ?????õĿհ????ݿ?outResult?????ã?ѭ???????и???ʼ
                    cbind(id=i,#?ϲ??У???ÿ????????ͳ??????
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],#??ȡ??????????HR
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],#??ȡ??????????HR??95%CI??ֵ
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],#??ȡ??????????HR??95%CI??ֵ
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#??ȡ??????????pֵ
    )
  }
}
#write.table(outResult,file="UniCoxSurvival.txt",sep="\t",row.names=F,quote=F) ##110 univaribles
outmRNAResult=data.frame() #??һ???հ????ݿ򣬺???forѭ????????
for(i in colnames(mRNA311.survival.train[,1:154])){ 
  tdcox <- coxph(Surv(Survival, Death) ~ mRNA311.survival.train[,i], data = mRNA311.survival.train)#??ʼ??һѭ??cox????
  tdcoxSummary = summary(tdcox) #summary??????tdcox?ܽᣬ??????????ȡ????
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #??ȡpֵ????????Ҫ?Ǻ?????ȡ????????gene??
  print(pvalue)
  print(i)
  if(pvalue< pFilter){ # ???????ǲ???Ҫ??ȡ???л????????ݣ?ֻ??Ҫ????????gene?Ľ?????????????????pvalue<0.05????ȡ????
    sigGenes=c(sigGenes,i)
    outmRNAResult=rbind(outmRNAResult,#?ϲ??У?ʵ?????Ƕ?ѭ???????ĺϲ???ǰ?????õĿհ????ݿ?outResult?????ã?ѭ???????и???ʼ
                    cbind(id=i,#?ϲ??У???ÿ????????ͳ??????
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],#??ȡ??????????HR
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],#??ȡ??????????HR??95%CI??ֵ
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],#??ȡ??????????HR??95%CI??ֵ
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#??ȡ??????????pֵ
    )
  }
}
#write.table(outmRNAResult,file="UniCoxSurvivalmRNA.txt",sep="\t",row.names=F,quote=F) ##110 univaribles

miRNA311.survival.train.uni<-miRNA311.survival.train[,outResult$id]
miRNA311.survival.train.uni<-cbind(miRNA311.survival.train$Survival,miRNA311.survival.train$Death,miRNA311.survival.train.uni)
colnames(miRNA311.survival.train.uni)[1]<-"Survival"
colnames(miRNA311.survival.train.uni)[2]<-"Death"
miRNA311.survival.train.uni.nozero<-miRNA311.survival.train.uni[which(miRNA311.survival.train.uni$Survival!=0),]

mRNA311.survival.train.uni<-mRNA311.survival.train[,outmRNAResult$id]
mRNA311.survival.train.uni<-cbind(mRNA311.survival.train$Survival,mRNA311.survival.train$Death,mRNA311.survival.train.uni)
colnames(mRNA311.survival.train.uni)[1]<-"Survival"
colnames(mRNA311.survival.train.uni)[2]<-"Death"

#scale2  <-  function(x, na.rm  =  FALSE)(x  -  mean(x, na.rm =  na.rm)) / sd(x, na.rm)

#####stepwise # use this not lasso cuz lasso 2 much variables
library(StepReg)

formula <-Surv(Survival,Death)~.
stepresult.AIC<-stepwiseCox(formula=formula,
            data=miRNA311.survival.train.uni.nozero,
            selection="bidirection",
           select="AIC", #could choose other parameter including AIC, AICc, SBC, IC(1), IC(3/2), HQ, HQc and Significant Levels(SL)
            sle=0.05,
            sls = 0.1
            )
#coef<-stepresult.AIC[[5]]$coef
#coef<-as.numeric(coef)
#names(coef)<-stepresult.AIC[[5]]$Variable


########Lasso
library(glmnet)
library(caret)
set.seed(1)

x <- as.matrix(miRNA311.survival.train.uni.nozero[,-c(1:2)])
x<-x[,colSums(is.infinite(x))==0] #rm inf columns
x<-apply(x, 2, function(x)as.numeric(x))
y<-data.matrix(Surv(miRNA311.survival.train.uni.nozero$Survival,miRNA311.survival.train.uni.nozero$Death))

xmRNA <- as.matrix(mRNA311.survival.train.uni[,-c(1:2)])
xmRNA<-xmRNA[,colSums(is.infinite(xmRNA))==0] #rm inf columns
xmRNA<-apply(xmRNA, 2, function(x)as.numeric(x))
ymRNA<-data.matrix(Surv(mRNA311.survival.train.uni$Survival,mRNA311.survival.train.uni$Death))


fit.miRNA <-glmnet(x,y,family = "cox",alpha = 1)
fit.mRNA<-glmnet(xmRNA,ymRNA,family = "cox",alpha = 1)
plot(fit.miRNA,label=T)
plot(fit.miRNA,xvar="lambda",label=T)
fitcv.miRNA <- cv.glmnet(x,y,family="cox", alpha=1,nfolds=10)
fitcv.mRNA <- cv.glmnet(xmRNA,ymRNA,family="cox", alpha=1,nfolds=10)

plot(fitcv.miRNA)
coef.miRNA<-coef(fitcv.miRNA,s="lambda.min")[which(coef(fitcv.miRNA, s="lambda.min")!=0),]
coef.mRNA<-coef(fitcv.mRNA,s="lambda.min")[which(coef(fitcv.mRNA, s="lambda.min")!=0),]
names(coef.miRNA)
names(coef.mRNA)
####barplot
#trellis.device(device="windows", height = 20, width = 40, color=TRUE)
par(mar=c(3,13,3,3))#????????
barplot(height = coef.miRNA[order(coef.miRNA)],  # ??ͼ???ݣ???ֵ????��??
        cex.names = 1.3,
        las=1, #??ת????????ǩ
        #family = 'Kai',  # ????????
       # axis(side=2,at = seq(-1, 2, .2)),
        # axisnames = T,
     #  border = '#ffffff', # ??��??ɫ
       # xlab = 'Ƶ??',  # X??????
       # ylab = '?Ա?',  # Y??????
      names=names(coef.miRNA[order(coef.miRNA)]),
        main = 'Coefficients',  # ??????
       horiz = TRUE,  # ?Ƿ?Ϊˮƽ????
      xlim = c(-0.5, 0.6), # X??ȡֵ??Χ
      col=c('orange','orange','orange','orange','orange','blue','blue','blue')
)

#####################calculate risk score##########################
###miRNA
miRNA.rs<-names(coef.miRNA)  #name of miRNA used to calculate risk score
miRNA.rs.train<-miRNA311.survival.train.uni.nozero[,miRNA.rs]
temp.train<-mapply(`*`,miRNA.rs.train, coef.miRNA)
riskscore.train<-rowSums(temp.train)
highrisk.miRNA<-which(riskscore.train>median(riskscore.train))
lowrisk.miRNA<-which(riskscore.train<=median(riskscore.train))

miRNA.rs.test<-miRNA311.survival.test[,miRNA.rs]
temp.test<-mapply(`*`,miRNA.rs.test, coef.miRNA)
riskscore.test<-rowSums(temp.test)

#highrisk.miRNA.test<-which(riskscore.test>median(riskscore.test))
#lowrisk.miRNA.test<-which(riskscore.test<=median(riskscore.test))

####mRNA
#mRNA.rs<-names(coef.mRNA)  #name of miRNA used to calculate risk score
#mRNA.rs.train<-mRNA311.survival.train.uni[,mRNA.rs]
#temp.train<-mapply(`*`,mRNA.rs.train, coef.mRNA)
#riskscore.train<-rowSums(temp.train)
#highrisk.mRNA<-which(riskscore.train>median(riskscore.train))
#lowrisk.mRNA<-which(riskscore.train<=median(riskscore.train))
#####forest plot
rownames(outResult)<-outResult$id
rsmiRNAuni<-outResult[miRNA.rs,]
colnames(rsmiRNAuni)<-c("Variable","HR_mean","Lower","Upper","P value")

rsmiRNAuni=rsmiRNAuni %>% mutate_at(c(2:5), as.numeric)
rsmiRNAuni<-data.frame(lapply(rsmiRNAuni, function(x)if(is.numeric(x))round(x,5) else x))

library(questionr)
write.csv(rsmiRNAuni,file="./rsmiRNAuni.csv")

######kaplan Meier curve train OS
miRNA311.survival.train.uni.nozero$risk<-"high"
miRNA311.survival.train.uni.nozero$risk[lowrisk.miRNA]<-"low"
mRNA311.survival.train.uni$risk<-"high"
mRNA311.survival.train.uni$risk[lowrisk.mRNA]<-"low"

fit <- survfit(Surv(Survival,Death) ~ risk,  # ???????????? 
               data = miRNA311.survival.train.uni.nozero) # ???ݼ?��Դ
OStrain<-as.numeric(surv_pvalue(fit)["pval"])
fit
summary(fit)
ggsurvplot(fit, # ?????????϶???
          # data = lung,  # ָ????��????��Դ
           conf.int = TRUE, # ??ʾ????????
           pval = TRUE, # ????Pֵ
           surv.median.line = "hv",  # ??????λ????ʱ????
           risk.table = TRUE, # ???ӷ??ձ?
           xlab = "Follow up time(Years)", # ָ??x????ǩ
           legend = c(0.9,0.9), # ָ??ͼ??λ??
           legend.title = "", # ????ͼ??????
           legend.labs = c("High-risk", "Low-risk"), # ָ??ͼ????????ǩ
           break.x.by = 2
           )  # ????x???̶ȼ???
##########RFS train
clinical<-read.csv("./Clinicalnature_Mstage_combine.csv")
survival.RFS<-clinical[,grep("RFS",colnames(clinical))]
survival.RFS<-survival.RFS[,-3]
colnames(survival.RFS)<-c("RFS","Recurrence")
survival.RFS$ID<-clinical$sampleID
rownames(survival.RFS)<-str_replace_all(clinical$sampleID,"-",".")
survival.RFS<-na.omit(survival.RFS)
RFS.train<-survival.RFS[rownames(miRNA311.survival.train.uni.nozero),]
RFS.train$risk<-miRNA311.survival.train.uni.nozero$risk
RFS.train$risk.mRNA<-mRNA311.survival.train.uni$risk
RFS.train<-na.omit(RFS.train)
RFS.train$RFS<-RFS.train$RFS/365
fit <- survfit(Surv(RFS,Recurrence) ~ risk,  # ???????????? risk or risk.mRNA
               data = RFS.train) # ???ݼ?��Դ
RFStrain<-as.numeric(surv_pvalue(fit)["pval"])

fit
summary(fit)
ggsurvplot(fit, # ?????????϶???
           # data = lung,  # ָ????��????��Դ
           conf.int = TRUE, # ??ʾ????????
           pval = TRUE, # ????Pֵ
           surv.median.line = "hv",  # ??????λ????ʱ????
           risk.table = TRUE, # ???ӷ??ձ?
           xlab = "Follow up time(Years)", # ָ??x????ǩ
           legend = c(0.9,0.9), # ָ??ͼ??λ??
           legend.title = "", # ????ͼ??????
           legend.labs = c("High-risk", "Low-risk"), # ָ??ͼ????????ǩ
           break.x.by = 2
)  # ????x???̶ȼ???

###########KM test OS
miRNA311.survival.test.nonzero<-miRNA311.survival.test[which(miRNA311.survival.test$Survival!=0),]
miRNA.rs.test<-miRNA311.survival.test.nonzero[,miRNA.rs]
temp.test<-mapply(`*`,miRNA.rs.test, coef.miRNA)
riskscore.test.miRNA<-rowSums(temp.test)

mRNA311.survival.test<-mRNA311.survival.test[which(mRNA311.survival.test$Survival!=0),]
mRNA.rs.test<-mRNA311.survival.test[,mRNA.rs]
temp.test<-mapply(`*`,mRNA.rs.test, coef.mRNA)
riskscore.test.mRNA<-rowSums(temp.test)


highrisk.test.miRNA<-which(riskscore.test.miRNA>median(riskscore.test.miRNA))
lowrisk.test.miRNA<-which(riskscore.test.miRNA<=median(riskscore.test.miRNA))

highrisk.test.mRNA<-which(riskscore.test.mRNA>median(riskscore.test.mRNA))
lowrisk.test.mRNA<-which(riskscore.test.mRNA<=median(riskscore.test.mRNA))


miRNA311.survival.test.nonzero$risk<-"high"
miRNA311.survival.test.nonzero$risk[lowrisk.test.miRNA]<-"low"

mRNA311.survival.test$risk.mRNA<-"high"
mRNA311.survival.test$risk.mRNA[lowrisk.test.mRNA]<-"low"

fit <- survfit(Surv(Survival,Death) ~ risk,  # ???????????? 
               data = miRNA311.survival.test.nonzero) # ???ݼ?��Դ
OStest<-as.numeric(surv_pvalue(fit)["pval"])

fit
summary(fit)
ggsurvplot(fit, # ?????????϶???
           # data = lung,  # ָ????��????��Դ
           conf.int = TRUE, # ??ʾ????????
           pval = TRUE, # ????Pֵ
           surv.median.line = "hv",  # ??????λ????ʱ????
           risk.table = TRUE, # ???ӷ??ձ?
           xlab = "Follow up time(Years)", # ָ??x????ǩ
           legend = c(0.9,0.9), # ָ??ͼ??λ??
           legend.title = "", # ????ͼ??????
           legend.labs = c("High-risk", "Low-risk"), # ָ??ͼ????????ǩ
           break.x.by = 2
)  # ????x???̶ȼ???
#########KM test RFS
RFS.test<-survival.RFS[rownames(miRNA311.survival.test.nonzero),]
RFS.test$risk<-miRNA311.survival.test.nonzero$risk
RFS.test<-na.omit(RFS.test)
RFS.test$RFS<-RFS.test$RFS/365
fit <- survfit(Surv(RFS,Recurrence) ~ risk,  # ???????????? 
               data = RFS.test) # ???ݼ?��Դ
RFStest<-as.numeric(surv_pvalue(fit)["pval"])

fit
summary(fit)
ggsurvplot(fit, # ?????????϶???
           # data = lung,  # ָ????��????��Դ
           conf.int = TRUE, # ??ʾ????????
           pval = TRUE, # ????Pֵ
           surv.median.line = "hv",  # ??????λ????ʱ????
           risk.table = TRUE, # ???ӷ??ձ?
           xlab = "Follow up time(Years)", # ָ??x????ǩ
           legend = c(0.9,0.9), # ָ??ͼ??λ??
           legend.title = "", # ????ͼ??????
           legend.labs = c("High-risk", "Low-risk"), # ָ??ͼ????????ǩ
           break.x.by = 2
)  # ????x???̶ȼ???
#####RFS.combine
miRNA.rs.RFS.combine<-miRNA311.survival[,miRNA.rs]
RFS.combine<-rbind(RFS.test,RFS.train[,-5])

miRNA.rs.RFS.combine<-miRNA.rs.RFS.combine[rownames(RFS.combine),]
temp.combine<-mapply(`*`,miRNA.rs.RFS.combine, coef.miRNA)
riskscore.RFS.combine<-rowSums(temp.combine)
fit <- survfit(Surv(RFS,Recurrence) ~ risk,  # ???????????? 
               data = RFS.combine) # ???ݼ?��Դ
RFScombine<-as.numeric(surv_pvalue(fit)["pval"])


fit
summary(fit)
ggsurvplot(fit, # ?????????϶???
           # data = lung,  # ָ????��????��Դ
           conf.int = TRUE, # ??ʾ????????
           pval = TRUE, # ????Pֵ
           surv.median.line = "hv",  # ??????λ????ʱ????
           risk.table = TRUE, # ???ӷ??ձ?
           xlab = "Follow up time(Years)", # ָ??x????ǩ
           legend = c(0.9,0.9), # ָ??ͼ??λ??
           legend.title = "", # ????ͼ??????
           legend.labs = c("High-risk", "Low-risk"), # ָ??ͼ????????ǩ
           break.x.by = 2
) 




p.df<-rbind(p.df,data.frame(t=t,OStrain=OStrain, RFStrain=RFStrain, OStest=OStest,RFStest=RFStest))

}
#########compare the miRNA level between nomal and tumor
miRNA.normal<-miRNA[,grep("11$",colnames(miRNA))]#tcga id end with 11 is normal
miRNA.tumor<-miRNA[,-grep("11$",colnames(miRNA))]
#miRNA.other<-miRNA[,grep("06",colnames(miRNA))]
miRNA.rs.ID<-str_replace_all(miRNA.rs,"\\.","\\-")#rs stands for risk score
miRNA.normal.rs<-miRNA.normal[miRNA.rs.ID,]
miRNA.normal.rs$miName<-rownames(miRNA.normal.rs)
miRNA.tumor.rs<-miRNA.tumor[miRNA.rs.ID,]
miRNA.tumor.rs$miName<-rownames(miRNA.tumor.rs)

boxplot.df.normal<-melt(miRNA.normal.rs,id.vars='miName')
boxplot.df.normal$Class<-"Normal"
boxplot.df.tumor<-melt(miRNA.tumor.rs,id.vars='miName')
boxplot.df.tumor$Class<-"Tumor"
boxplot.df<-rbind(boxplot.df.normal,boxplot.df.tumor)
boxplot.df$logvalue<-log2(boxplot.df$value)
df_p_val <- boxplot.df %>% 
  group_by(miName) %>% 
  wilcox_test(formula = value ~ Class) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='miRNA')
#p1 <- ggplot(boxplot.df,aes(x=miName,y=logvalue,fill=Class))+
#  geom_boxplot(width=0.6,alpha=0.8)
#win.graph(width=4.875, height=2.5,pointsize=8)
ggplot()+
  geom_boxplot(data = boxplot.df,mapping = aes(x=miName,y=value,fill=Class),width=0.5)+
  scale_fill_manual(values = c('#E21C21','#3A7CB5'))+
 stat_pvalue_manual(
    df_p_val, x = "miName", y.position = df_p_val$y.position+2,
    label = "p.signif",label.size=3,
    position = position_dodge(0.8)
  )+
  labs(y='miRNA Expression (Log10)',x='')+
  theme_test()+
  theme(axis.text = element_text(color = 'black',size = 15),#????????????С
        plot.caption = element_markdown(face = 'bold'),
        legend.position = c(0.85,0.95),legend.text = element_text (size = 15),#ͼ????????С
        legend.title = element_text(size=15),#ͼ???????ִ?С
        axis.title.y=element_text(size=15),#y??label??С
        legend.direction = 'horizontal')+#ͼ??????
  theme (axis.text.x = element_text (angle=45, hjust=1, vjust=1))
#####heatmap visualize the difference expression between highrisk and lowrisk
library(pheatmap)
matrix.miRNA<-miRNA311.survival[,miRNA.rs]
matrix.miRNA<-as.data.frame(t(matrix.miRNA))
temp.all<-mapply(`*`,as.data.frame(t(matrix.miRNA)), coef.miRNA)
riskscore.all<-rowSums(temp.all)
highrisk<-which(riskscore.all>median(riskscore.all))
lowrisk<-which(riskscore.all<=median(riskscore.all))

annotation<-data.frame(Sample=rownames(miRNA311.survival),Risk="High risk")
annotation$Risk[lowrisk]<-"Low risk"
rownames(annotation)<-annotation$Sample
group1=annotation[colnames(matrix.miRNA),-1,drop=FALSE]
group.risk<-annotation[order(group1$Risk),-1,drop=FALSE]
#trellis.device(device="windows", height = 6, width = 25, color=TRUE)
matrix.miRNA.df<-matrix.miRNA[,order(group1$Risk)]
pheatmap(matrix.miRNA.df,
         cluster_rows = FALSE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         scale="row",
         cluster_cols = FALSE,
         fontsize_row = 10,
         fontsize_col = 10,
         annotation_col = group.risk,
         #color = colorRampPalette(c("navy", "white", "red"))(100),
         color = colorRampPalette(c("red", "white", "blue4"))(100),#????ɫ
         angle_col = 45, #?޸ĺ???????????б??
        # filename = 'cor.fpkm2.png',
)
#####spearman correlation
scatterplot.data<-t(matrix.miRNA)
scatterplot.data<-as.data.frame(scatterplot.data)
scatterplot.data$riskscore<-riskscore.all
scatterplot.data$Risk<-annotation$Risk
correls(riskscore.all, t(matrix.miRNA), type = "spearman")
library(patchwork)

h<- ggplot(scatterplot.data, aes(scatterplot.data[,1],scatterplot.data[,9]))+
    geom_point()+geom_smooth(method="loess",color="blue",fill="lightgray")+
  #stat_cor(method="spearman",label.y=4.5,label.x=4)
    stat_cor(method="spearman",label.y.npc = "top",label.x.npc ="left")+
    xlab("hsa-miR-106b-5p")+ylab("Riskscore")
  
b<- ggplot(scatterplot.data, aes(scatterplot.data[,2],scatterplot.data[,9]))+
  geom_point()+geom_smooth(method="loess",color="blue",fill="lightgray")+
  stat_cor(method="spearman",label.y.npc = "top",label.x.npc ="left")+
 xlab("hsa-miR-3142")+ylab("Riskscore")

c<- ggplot(scatterplot.data, aes(scatterplot.data[,3],scatterplot.data[,9]))+
  geom_point()+geom_smooth(method="loess",color="blue",fill="lightgray")+
  stat_cor(method="spearman",label.y.npc = "top",label.x.npc ="left")+
  xlab("hsa-miR-3176")+ylab("Riskscore")

d<- ggplot(scatterplot.data, aes(scatterplot.data[,4],scatterplot.data[,9]))+
  geom_point()+geom_smooth(method="loess",color="blue",fill="lightgray")+
  stat_cor(method="spearman",label.y.npc = "top",label.x.npc ="left")+
  xlab("hsa-miR-5707")+ylab("Riskscore")

e<- ggplot(scatterplot.data, aes(scatterplot.data[,5],scatterplot.data[,9]))+
  geom_point()+geom_smooth(method="loess",color="blue",fill="lightgray")+
  stat_cor(method="spearman",label.y.npc = "top",label.x.npc ="left")+
  xlab("hsa-miR-6865-5p")+ylab("Riskscore")

f<- ggplot(scatterplot.data, aes(scatterplot.data[,6],scatterplot.data[,9]))+
  geom_point()+geom_smooth(method="loess",color="blue",fill="lightgray")+
  stat_cor(method="spearman",label.y.npc = "top",label.x.npc ="left")+
  xlab("hsa-miR-7109-5p")+ylab("Riskscore")

g<- ggplot(scatterplot.data, aes(scatterplot.data[,7],scatterplot.data[,9]))+
  geom_point()+geom_smooth(method="loess",color="blue",fill="lightgray")+
  stat_cor(method="spearman",label.y.npc = "top",label.x.npc ="middle")+
  xlab("hsa-miR-769-3p")+ylab("Riskscore")

a<- ggplot(scatterplot.data, aes(scatterplot.data[,8],scatterplot.data[,9]))+
  geom_point()+geom_smooth(method="loess",color="blue",fill="lightgray")+
  stat_cor(method="spearman",label.y.npc = "top",label.x.npc ="left")+
  xlab("hsa-miR-96-5p")+ylab("Riskscore")

a+b+c+d+e+f+g+h
#####################Figure2
par(mar=c(5,5,3,3))

riskscore.train.plot<-data.frame(riskscore=riskscore.train,group=1)
riskscore.train.plot$ID<-rownames(miRNA311.survival.train.uni.nozero)
riskscore.train.plot$Death<-miRNA311.survival.train.uni.nozero$Death
riskscore.train.plot$Death[which(riskscore.train.plot$Death==1)]<-2
riskscore.train.plot$Death[which(riskscore.train.plot$Death==0)]<-1

riskscore.train.plot$Survival<-miRNA311.survival.train.uni.nozero$Survival

riskscore.train.plot<-riskscore.train.plot[order(riskscore.train),]
riskscore.train.plot$group[which(riskscore.train.plot$riskscore>median(riskscore.train))]<-"2"
split.screen(c(2, 1))
screen(1)
plot(riskscore.train.plot$riskscore,col=riskscore.train.plot$group,ylab="Riskscore",xlab="Patient(increasing risk score)")
abline(v=295,col = 'coral2', lwd = 2,lty=2)
screen(2)
plot(riskscore.train.plot$Survival,col=riskscore.train.plot$Death)
wilcox.test(riskscore.train.plot$riskscore[which(riskscore.train.plot$Death==1)],riskscore.train.plot$riskscore[which(riskscore.train.plot$Death==2)])
####test
riskscore.test.plot<-data.frame(riskscore=riskscore.test.miRNA,group=1)
riskscore.test.plot$ID<-rownames(miRNA311.survival.test.nonzero)
riskscore.test.plot$Death<-miRNA311.survival.test.nonzero$Death
riskscore.test.plot$Death[which(riskscore.test.plot$Death==1)]<-2
riskscore.test.plot$Death[which(riskscore.test.plot$Death==0)]<-1
#riskscore.test.plot$Deathori<-miRNA311.survival.test.nonzero$Death

riskscore.test.plot$Survival<-miRNA311.survival.test.nonzero$Survival

riskscore.test.plot<-riskscore.test.plot[order(riskscore.test),]
riskscore.test.plot$group[which(riskscore.test.plot$riskscore>median(riskscore.test))]<-"2"
split.screen(c(2, 1))
screen(1)
plot(riskscore.test.plot$riskscore,col=riskscore.test.plot$group,ylab="Riskscore",xlab="Patient(increasing risk score)")
abline(v=295,col = 'coral2', lwd = 2,lty=2)
screen(2)
plot(riskscore.test.plot$Survival,col=riskscore.test.plot$Death)
source("H:/COXEN/breast/coxen-library-2016-10-25.R")
riskscore.all.plot<-rbind(riskscore.train.plot,riskscore.test.plot)
wilcox.test(riskscore.all.plot$riskscore[which(riskscore.all.plot$Death==1)],riskscore.all.plot$riskscore[which(riskscore.all.plot$Death==2)],alternative ="two.sided",correct = FALSE)
########RFS.combine
RFS.combine$riskscore<-riskscore.RFS.combine
RFS.combine$group<-"1"
RFS.combine$group[which(RFS.combine$riskscore>median(RFS.combine$riskscore))]<-"2"
RFS.combine$Recurrence[which(RFS.combine$Recurrence==1)]<-"2"
RFS.combine$Recurrence[which(RFS.combine$Recurrence==0)]<-"1"
RFS.combine<-RFS.combine[order(RFS.combine$riskscore),]

split.screen(c(2, 1))
screen(1)
plot(RFS.combine$riskscore,col=RFS.combine$group,ylab="Riskscore",xlab="Patient(increasing risk score)")
abline(v=295,col = 'coral2', lwd = 2,lty=2)
screen(2)
plot(RFS.combine$RFS,col=RFS.combine$Recurrence)

##########################
#plot.strip( a=riskscore.test.plot$riskscore, b=riskscore.test.plot$Deathori)
#boxplot(riskscore.test.plot$riskscore~riskscore.test.plot$Death)
mycomparison=list(c("1","2"))
ggviolin(riskscore.all.plot, x = "Death", y = "riskscore", fill = "Death",legend= "none",palette =c("#00AFBB","#FC4E07"),
         add = "boxplot",xlab="Living status",ylab="Riskscore",add.params = list(fill="white"))+
         stat_compare_means(comparisons = mycomparison,
                           # method = "wilcox.test",
                            #label="p.signif",  #???Ǻ???ʾpֵ
                            #aes(label = "p.format"),#??ʵ?ʵ?pֵ
                            label.y =2)+
        scale_x_discrete(labels = c("1" = "Live","2" = "Dead"))

####train.clinical to do multivariate analysis
library(stringr)
train.clinical<-RFS.train
train.clinical$OS<-miRNA311.survival.train.uni.nozero[rownames(train.clinical),]$Survival
train.clinical$Death<-miRNA311.survival.train.uni.nozero[rownames(train.clinical),]$Death
rownames(clinical)<-str_replace_all(clinical$sampleID,"-",".")
train.clinical$Stage<-clinical[rownames(train.clinical),]$pathologic_stage
train.clinical$M<-clinical[rownames(train.clinical),]$pathologic_M
train.clinical$T<-clinical[rownames(train.clinical),]$pathologic_T
train.clinical$N<-clinical[rownames(train.clinical),]$pathologic_N
train.clinical$Age<-clinical[rownames(train.clinical),]$age_at_initial_pathologic_diagnosis
train.clinical$PAM50<-clinical[rownames(train.clinical),]$PAM50Call_RNAseq
train.clinical$Lymphnodes<-clinical[rownames(train.clinical),]$number_of_lymphnodes_positive_by_he
####test
test.clinical<-RFS.test
test.clinical$OS<-miRNA311.survival.test.nonzero[rownames(test.clinical),]$Survival
test.clinical$Death<-miRNA311.survival.test.nonzero[rownames(test.clinical),]$Death
rownames(clinical)<-str_replace_all(clinical$sampleID,"-",".")
test.clinical$Stage<-clinical[rownames(test.clinical),]$pathologic_stage
test.clinical$M<-clinical[rownames(test.clinical),]$pathologic_M
test.clinical$T<-clinical[rownames(test.clinical),]$pathologic_T
test.clinical$N<-clinical[rownames(test.clinical),]$pathologic_N
test.clinical$Age<-clinical[rownames(test.clinical),]$age_at_initial_pathologic_diagnosis
test.clinical$PAM50<-clinical[rownames(test.clinical),]$PAM50Call_RNAseq
test.clinical$Lymphnodes<-clinical[rownames(test.clinical),]$number_of_lymphnodes_positive_by_he
rownames(riskscore.test.plot)<-riskscore.test.plot$ID
test.clinical$riskscore<-riskscore.test.plot[rownames(test.clinical),]$riskscore
test.clinical.multi<-test.clinical

test.clinical.multi<-within(test.clinical.multi,rm(RFS,Recurrence))
#####
library(survival)
rownames(riskscore.train.plot)<-riskscore.train.plot$ID
train.clinical$riskscore<-riskscore.train.plot[rownames(train.clinical),]$riskscore

train.clinical.multi<-train.clinical
train.clinical.multi<-na.omit(train.clinical)

train.clinical.multi<-within(train.clinical.multi,rm(Recurrence,RFS,risk.mRNA))

train.clinical.multi$Stage[grep ("Stage I$|Stage IA|Stage IB",train.clinical.multi$Stage)]<-1
train.clinical.multi$Stage[grep ("Stage II$|Stage IIA|Stage IIB",train.clinical.multi$Stage)]<-2
train.clinical.multi$Stage[grep ("Stage III$|Stage IIIA|Stage IIIB|Stage IIIC",train.clinical.multi$Stage)]<-3
train.clinical.multi$Stage[grep ("Stage IV$",train.clinical.multi$Stage)]<-4
train.clinical.multi$Stage[grep ("Stage X$",train.clinical.multi$Stage)]<-5
#train.clinical.multi<-train.clinical.multi[-which(train.clinical.multi$Stage=="[Discrepancy]"),]
train.clinical.multi$Stage[which(train.clinical.multi$Stage=="")]<-NA

train.clinical.multi$M[grep("cM0",train.clinical.multi$M)]<-"M0"
train.clinical.multi$N[grep("N0",train.clinical.multi$N)]<-1
train.clinical.multi$N[grep("N1",train.clinical.multi$N)]<-1
train.clinical.multi$N[grep("N2",train.clinical.multi$N)]<-2
train.clinical.multi$N[grep("N3",train.clinical.multi$N)]<-2

train.clinical.multi$T[grep("T1",train.clinical.multi$T)]<-1
train.clinical.multi$T[grep("T2",train.clinical.multi$T)]<-1
train.clinical.multi$T[grep("T3",train.clinical.multi$T)]<-2
train.clinical.multi$T[grep("T4",train.clinical.multi$T)]<-2
train.clinical.multi$Stage<-as.factor(train.clinical.multi$Stage)

#train.clinical.multi<-train.clinical.multi[,-c(1,9)]
fit <- coxph( Surv(OS, Death) ~ riskscore+Stage+Age+T+N+M, data=train.clinical.multi)

#temp<-train.clinical.multi[!is.na(train.clinical.multi$Stage),]
####combined clinical to do multivariate analysis
clinical.multi<-rbind(within(train.clinical,rm(risk.mRNA)),test.clinical)
clinical.multi$Stage[grep ("Stage I$|Stage IA|Stage IB",clinical.multi$Stage)]<-1
clinical.multi$Stage[grep ("Stage II$|Stage IIA|Stage IIB",clinical.multi$Stage)]<-2
clinical.multi$Stage[grep ("Stage III$|Stage IIIA|Stage IIIB|Stage IIIC",clinical.multi$Stage)]<-3
clinical.multi$Stage[grep ("Stage IV",clinical.multi$Stage)]<-4
clinical.multi$Stage[which(clinical.multi$Stage=="")]<-NA
clinical.multi$Stage<-as.numeric(clinical.multi$Stage)

clinical.multi$M[grep("cM0",clinical.multi$M)]<-"NA"
clinical.multi$M[grep("M1",clinical.multi$M)]<-"1"
clinical.multi$M[grep("M0",clinical.multi$M)]<-"0"
clinical.multi$M[grep("MX",clinical.multi$M)]<-"NA"

clinical.multi$M<-as.numeric(clinical.multi$M)

clinical.multi$N[grep("N0",clinical.multi$N)]<-0
clinical.multi$N[grep("N1",clinical.multi$N)]<-1
clinical.multi$N[grep("N2",clinical.multi$N)]<-2
clinical.multi$N[grep("N3",clinical.multi$N)]<-3
clinical.multi$N<-as.numeric(clinical.multi$N)

clinical.multi$T[grep("T1|TX",clinical.multi$T)]<-1
clinical.multi$T[grep("T2",clinical.multi$T)]<-1
clinical.multi$T[grep("T3",clinical.multi$T)]<-2
clinical.multi$T[grep("T4",clinical.multi$T)]<-2
clinical.multi$T<-as.numeric(clinical.multi$T)

clinical.multi$PAM50[grep("Normal",clinical.multi$PAM50)]<-1
clinical.multi$PAM50[grep("LumB",clinical.multi$PAM50)]<-2
clinical.multi$PAM50[grep("LumA",clinical.multi$PAM50)]<-3
clinical.multi$PAM50[grep("Her2",clinical.multi$PAM50)]<-4
clinical.multi$PAM50[grep("Basal",clinical.multi$PAM50)]<-5
clinical.multi$PAM50<-as.factor(clinical.multi$PAM50)

fit <- coxph( Surv(OS, Death) ~Stage, data=clinical.multi)
fit <- coxph( Surv(OS, Death) ~ riskscore+Stage+Age+M, data=clinical.multi)
######forestplot
df.forestplot<-read.csv("./multiforestplot.csv")
colnames(df.forestplot)[1]<-"HR_mean"
library(forestplot)

forestplot(labeltext=as.matrix(df.forestplot[,c(5,6,4)]),
           mean=df.forestplot$HR_mean,
           lower=df.forestplot$Lower,
           upper=df.forestplot$Upper,
           #is.summary=c(T,T,T,F,F,T,F,F,T,T,T,T,T,T,F,F,T,F,F,T,F,F,T,T,T,T,T),
           graph.pos=2, 
           zero=1, 
           graphwidth=unit(40,"mm"),
           lwd.ci = 1.5,# ?????????ߵĿ???
           ci.vertices.height = 0.02, # # ??????ĩ?˵ĳ???
           lineheight="auto", #?ߵĸ߶?eg:unit(20,"cm")
           line.margin=0.1,
           boxsize=0.05,  ##???????е?Բ?ĵ???С
           # xlog=TRUE,#x????????ȡ????
           #   xlim=c(0,50),
           lty.ci = 7,# ?????????ߵ?????
           xticks.digits = 1,
           xticks=c(1,5,10,20,50),
           fn.ci_norm = fpDrawCircleCI,#????????ʾ??ʽ
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.6), xlab = gpar(cex = 0.8), cex = 0.9), #?ı???С
           col= fpColors(line = "#CC79A7", #?????????ߵ?????
                         box="#D55E00")#?????????ߵĿ???
)


#####nomogram
dd <- datadist(clinical.multi)
options(datadist = "dd")
cph<-cph(Surv(OS, Death)~Stage+Age+riskscore,data = clinical.multi,x = TRUE, y = TRUE, surv = TRUE)

survival <- Survival(cph)
survival1 <- function(x) survival(1, x)
survival3 <- function(x) survival(3, x)
survival5 <- function(x) survival(5, x)
survival10 <- function(x) survival(10, x)

nom <- nomogram(cph, fun = list(survival1, survival3,survival5), fun.at = c(0.05, seq(0.1,
      0.9, by = 0.1), 0.95), funlabel = c("1 year survival", "3 year survival","5 year survival"))
plot(nom)
###ROC ????
#???ɷ?????Ҫ????????ʽ?????????ݷֱ?Ϊ????ʱ?䡢????״̬??CoxԤ??ֵ
pred_f_training.nomo<-predict(cph,clinical.multi,type="lp")#!!!type="lp",????û??
pred_f_training.stage<-as.numeric(clinical.multi$Stage)
pred_f_training.riskscore<-as.numeric(clinical.multi$riskscore)

data_table<-data.frame(time=clinical.multi[,"OS"],status=clinical.multi[,"Death"],score.nomo=pred_f_training.nomo,score.stage=pred_f_training.stage,score.riskscore=pred_f_training.riskscore)
#?ⲿ??????Ԥ?ⲻ??1?ꡢ3?ꡢ5???Ļ?????Ҫ????times?????е???ֵ
library(timeROC)
time_roc_stage <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score.stage,
  cause = 1,
  weighting="marginal",
  times = c(1, 3, 5),
  ROC = TRUE,
  iid = TRUE
)
time_roc_riskscore <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score.riskscore,
  cause = 1,
  weighting="marginal",
  times = c(1, 3, 5),
  ROC = TRUE,
  iid = TRUE
)
time_roc_nomo <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score.nomo,
  cause = 1,
  weighting="marginal",
  times = c(1, 3, 5),
  ROC = TRUE,
  iid = TRUE
)
AUC.riskscore<-time_roc_riskscore$AUC
AUC.Stage<-time_roc_stage$AUC
AUC.nomo<-time_roc_nomo$AUC
time_ROC_nomo <- data.frame(
  TP_1year = time_roc_nomo$TP[, 1],
  FP_1year = time_roc_nomo$FP[, 1],
  TP_3year = time_roc_nomo$TP[, 2],
  FP_3year = time_roc_nomo$FP[, 2],
  TP_5year = time_roc_nomo$TP[, 3],
  FP_5year = time_roc_nomo$FP[, 3]
)
time_ROC_stage <- data.frame(
  TP_1year = time_roc_stage$TP[, 1],
  FP_1year = time_roc_stage$FP[, 1],
  TP_3year = time_roc_stage$TP[, 2],
  FP_3year = time_roc_stage$FP[, 2],
  TP_5year = time_roc_stage$TP[, 3],
  FP_5year = time_roc_stage$FP[, 3]
)
time_ROC_riskscore <- data.frame(
  TP_1year = time_roc_riskscore$TP[, 1],
  FP_1year = time_roc_riskscore$FP[, 1],
  TP_3year = time_roc_riskscore$TP[, 2],
  FP_3year = time_roc_riskscore$FP[, 2],
  TP_5year = time_roc_riskscore$TP[, 3],
  FP_5year = time_roc_riskscore$FP[, 3]
)
#colnames(time_ROC_riskscore)<-paste0(colnames(time_ROC_riskscore),"_","riskscore")
#colnames(time_ROC_Stage)<-paste0(colnames(time_ROC_Stage),"_","Stage")
#colnames(time_ROC_nomo)<-paste0(colnames(time_ROC_nomo),"_","nomo")
#time_roc_all<-list("nomo"=time_ROC_nomo,"riskscore"=time_ROC_riskscore,"Stage"=time_ROC_Stage)


#????1??3??5????????Ԥ????ROCͼ

ggplot(data = time_ROC_riskscore) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text", x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 years = ", 
                          sprintf("%.3f", time_roc_riskscore$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text", x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", 
                          sprintf("%.3f", time_roc_riskscore$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text", x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", 
                          sprintf("%.3f", time_roc_riskscore$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", 
                                margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", 
                                margin = margin(c(0, 15, 0, 0)))
  )


####calibration curve
units(clinical.multi$OS)<-"Year" ###define the unit of time, Year or Month or Day
cph1<-cph(Surv(OS, Death)~Stage+Age+riskscore,data = clinical.multi,x = TRUE, y = TRUE, surv = TRUE,time.inc=3)

cal1 <- calibrate(cph1, cmethod = "KM", method = "boot", u=3,m = 70, B = 1000)

plot(cal1, lwd = 2, lty = 0,
     errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     xlab = "Nomogram-Predicted Probability of 3-Year OS", 
     ylab = "Actual 3-Year DFS (proportion)",
     col = c(rgb(192, 98, 83, maxColorValue = 255)), 
     subtitles = FALSE, xlim = c(0.7,1), ylim = c(0.7, 1),
     main = "Calibrate plot")

                                                                             


################subtype compare riskscore##################
clinical.multi<-rbind(within(train.clinical,rm(risk.mRNA)),test.clinical)
clinical.nature.pam50<-clinical[rownames(clinical.multi),]
clinical.multi$pamNature<-clinical.nature.pam50$PAM50_mRNA_nature2012
clinical.multi$Her2<-clinical.nature.pam50$HER2_Final_Status_nature2012
clinical.multi$ER<-clinical.nature.pam50$ER_Status_nature2012
clinical.multi$PR<-clinical.nature.pam50$PR_Status_nature2012
clinical.multi$TNBC<-"Other"
clinical.multi$TNBC[which(clinical.multi$PR=="Negative" & clinical.multi$ER=="Negative" & clinical.multi$Her2=="Negative")]<-"TNBC"
#clinical.multi.nona<-clinical.multi[!is.na(clinical.multi$pamNature),]
clinical.multi.nona<-clinical.multi[complete.cases(clinical.multi[,c("PR","ER","Her2")]),]

#clinical.multi.pam50<-clinical.multi[which(clinical.multi$PAM50 !="Normal"),] 

boxplot(riskscore~TNBC,data=clinical.multi.nona)

# 获取 riskscore 的最小值和最大值
min_value <- min(clinical.multi.nona$riskscore)
max_value <- max(clinical.multi.nona$riskscore)

# 计算 y 轴间隔
interval <- 0.2

# 将 y 轴最小值和最大值调整为最接近的间隔值
min_value <- floor(min_value / interval) * interval
#max_value <- ceiling(max_value / interval) * interval
max_value <- 0.4

# 创建 ggviolin 图并添加统计检验
ggviolin(clinical.multi.nona, 
         x = "TNBC", 
         y = "riskscore", 
         fill = "TNBC",
         legend = "none",
         palette = "jco",
         add = "boxplot",
         xlab = "Living status",
         ylab = "Riskscore",
         add.params = list(fill = "white")) +
  stat_compare_means(comparisons = mycomparison,
                     method = "wilcox.test",
                     label.y = 0.4) +
  scale_y_continuous(breaks = seq(min_value, max_value, by = interval))
######validation
############GSE19783 
miRNA.ILMN.ID<-c("ILMN_3167226","ILMN_3167865","ILMN_3168022")#no3142,no 3176,no 5707,no 6865,no 7109,
miRNA.ID<-c("hsa-miR-96" ,"hsa-miR-106b","hsa-miR-769-3p")#no3142,no 3176,no 5707,no 6865,no 7109,

matrix.19783<-read.csv("./GSE19783.matrix.csv")
clinical.19783<-read.csv("./GSE19783.clinical.csv")
#clinical.19783<-na.omit(clinical.19783)
rownames(clinical.19783)<-clinical.19783$ID
rownames(matrix.19783)<-matrix.19783$ID
#matrix.19783<-matrix.19783[,rownames(clinical.19783)]
sel.matrix<-matrix.19783[miRNA.ID,]
sel.matrix<-t(sel.matrix)
sel.matrix<-as.data.frame(sel.matrix)
sel.matrix<-sel.matrix[-1,]
sel.matrix<-sel.matrix %>% mutate_at(c(1:3), as.numeric)

sel.coef.miRNA<-coef.miRNA[c(1:2,8)]
product<-mapply(`*`,sel.matrix, sel.coef.miRNA)
riskscore.19783<-rowSums(product)
highrisk.19783<-which(riskscore.19783>median(riskscore.19783))
lowrisk.19783<-which(riskscore.19783<=median(riskscore.19783))

rownames(clinical.19783)<-clinical.19783$ID
clinical.19783$risk<-"high"
clinical.19783$risk[lowrisk.19783]<-"low"
clinical.19783$DFS<-clinical.19783$DFS/12
clinical.19783$Status[which(clinical.19783$Status=="Dead")]<-"1"
clinical.19783$Status[which(clinical.19783$Status=="Alive")]<-"0"

clinical.19783<-na.omit(clinical.19783)
clinical.19783$Status<-as.numeric(clinical.19783$Status)

clinical.19783.sub<-clinical.19783[which(clinical.19783$Subtypes !=""),]
clinical.19783.her2<-clinical.19783[which(clinical.19783$X.Sample_characteristics_ch1.2 !=""),]
fit <- survfit(Surv(DFS,Status) ~ risk,  
               data = clinical.19783.her2)
dfs.19783<-as.numeric(surv_pvalue(fit)["pval"])

fit
summary(fit)
ggsurvplot(fit,
           conf.int = TRUE, 
           pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, 
           xlab = "Follow up time(Years)",
           legend = c(0.9,0.9), 
           legend.title = "", 
           legend.labs = c("High-risk", "Low-risk"),
           break.x.by = 2
)  
########GSE22220

miRNA.ILMN.ID<-c("ILMN_3167226","ILMN_3167865","ILMN_3168022")#no3142,no 3176,no 5707,no 6865,no 7109,

matrix.22220<-read.csv("./GSE22220_matrix.csv")
clinical.22220<-read.csv("./GSE22220_clinical.csv")
rownames(matrix.22220)<-matrix.22220$ID
sel.matrix<-matrix.22220[miRNA.ILMN.ID,]
sel.matrix<-t(sel.matrix)
sel.matrix<-as.data.frame(sel.matrix)
sel.matrix<-sel.matrix[-1,]
sel.matrix<-sel.matrix %>% mutate_at(c(1:3), as.numeric)

sel.coef.miRNA<-coef.miRNA[c(1:2,8)]
product<-mapply(`*`,sel.matrix, sel.coef.miRNA)
riskscore.22220<-rowSums(product)
highrisk.22220<-which(riskscore.22220>median(riskscore.22220))
lowrisk.22220<-which(riskscore.22220<=median(riskscore.22220))

rownames(clinical.22220)<-clinical.22220$ID
clinical.22220$risk<-"high"
clinical.22220$risk[lowrisk.22220]<-"low"
clinical.22220$DFS<-clinical.22220$DFS/12
#clinical.22220$Status<-as.factor(clinical.22220$Status)
clinical.22220<-na.omit(clinical.22220)
fit <- survfit(Surv(DRFS,Status) ~ risk,  
               data = clinical.22220)
dfs.22220<-as.numeric(surv_pvalue(fit)["pval"])

fit
summary(fit)
ggsurvplot(fit,
           conf.int = TRUE, 
           pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, 
           xlab = "Follow up time(Years)",
           legend = c(0.9,0.9), 
           legend.title = "", 
           legend.labs = c("High-risk", "Low-risk"),
           break.x.by = 2
)  

#########waterfall mutation between hign and low risk groups
library(maftools)
highrisk.tcga.ID<-clinical.multi$ID[which(clinical.multi$risk=="high")]
lowrisk.tcga.ID<-clinical.multi$ID[which(clinical.multi$risk=="low")]
annotate333<-data.frame(ID=c(highrisk.tcga.ID,lowrisk.tcga.ID),class=c(rep("high",166),rep("low",167)))
rownames(annotate333)<-c(highrisk.tcga.ID,lowrisk.tcga.ID)
tcga.mutation<-read.csv("./tcga_mutation.csv")
tcga.mutation$Sample_ID12<-str_replace_all(tcga.mutation$Sample_ID,".$","")
tcga.mutation.high<-tcga.mutation[tcga.mutation$Sample_ID12 %in% highrisk.tcga.ID,]
tcga.mutation.low<-tcga.mutation[tcga.mutation$Sample_ID12 %in% lowrisk.tcga.ID,]
tcga.mutation.low<-tcga.mutation.low[,-12]
tcga.mutation.high<-tcga.mutation.high[,-12]
write.table(tcga.mutation.high, 'tcga.mutation.high.tsv', sep ="\t",col.names = NA)
write.table(tcga.mutation.low, 'tcga.mutation.low.tsv', sep ="\t",col.names = NA)

options(stringsAsFactors = F) 
library(data.table)
tmp=fread('tcga.mutation.low.tsv')
tmp<-tmp[,-1]
head(tmp)   
colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tail(sort(table(tmp$Variant_Classification )))
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
table(tmp$Variant_Type )
tcga.brca = read.maf(maf = tmp,
                     vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))

oncoplot(maf = tcga.brca, top = 20) # 高频突变的前10个基因
gene_freq<-table(tcga.brca@data$Hugo_Symbol)
head(gene_freq[order(-gene_freq)])



#######differential expressed genes

library(data.table)
#?data.table
#R中的data.table包提供了一个data.frame的高级版本，让你的程序做数据整型的运算速度大大的增加。
a <- fread("./UCSC-TCGA-BRCA.htseq_counts.tsv.gz",sep = '\t',header = T)
#读取表达矩阵文件
a <- as.data.frame(a)
#转换为数据框，读取后，转换前的a的类型是"data.table" "data.frame"
a[1:4,1:4]
rownames(a) <- a[,1]
a <- a[,-1]
#去除第一列，行名修改。
genes <- row.names(a)
genes[1:10]
##接下来还原counts数，网站使用了log2(count+1)进行counts数转换，接下来进行还原
a <- 2^a-1
colnames(a)<-str_replace_all(colnames(a),".$","")
#colnames(a)<-str_replace_all(colnames(a),"\\-","\\.")
exprSet <- a[,intersect(rownames(annotate333),colnames(a))]
annnote332<-annotate333[intersect(rownames(annotate333),colnames(a)),]
table(apply(exprSet, 1, function(x){sum(x==0)<160}))
tmp <- apply(exprSet, 1, function(x){sum(x==0)<160})
exprSet <- ceiling(exprSet)

#对exprSet进行按行执行函数function：
#如果表达矩阵是0则取1，在18个样本中，如果这个基因的表达矩阵超过160个样本都是空值，则舍弃这个基因
#tmp是一个logical向量
exprSet <- exprSet[tmp,]
library(DESeq2)
#构建一个病例号和肿瘤分类的对应关系
colData <- data.frame(row.names = colnames(exprSet),group_list= annnote332$class)
#colData:是患者ID号和组织类型的对应关系

#构建DESeq()函数要求的表达式
dds <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData,design = ~ group_list)
#countData = exprSet,指出DESeq的表达矩阵
#colData = colData,知名每个表达矩阵的分类，比如实验组&对照组，正常组织&癌症组织
#question
#design= ~group_list，因子型，指出不同组的区别，是有顺序的。
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~group_list)
dds <- DESeq(dds)
#进行差异基因分析
resultsNames(dds)
res <-  results(dds, contrast=c("group_list","low","high"))
#用group_list来做引导文件，用tumor来比较normal组织

resOrdered <- res[order(res$padj),]
#把res差异分析文件通过padj来排序
head(resOrdered)
resOrdered=as.data.frame(resOrdered)
#把resOrdered变成数据框，
write.table(resOrdered, file="BRCAhigh-low.xls", sep="\t",quote=F)


## 我们使用|logFC| > 0.5，padj < 0.05（矫正后P值）
#LogFC =1
padj = 0.05
## 筛选出所有差异基因的结果
All_diffSig <- resOrdered[(resOrdered$padj < padj & (resOrdered$log2FoldChange>1 | resOrdered$log2FoldChange < (-1))),]
#---------------------
dim(All_diffSig)
[1] 2582 6

## 我们发现竟然没有差异基因，这是应该我这边的数据是随机的结果，如果你的数据有这样的问题，你需要在仔细检查一下哦。
## 我们为了下面的操作正常进行，我们选用的P值（未矫正）进行筛选。
write.csv(All_diffSig, "all.diffsig.csv")  ##输出差异基因数据集
#筛选上调和下调基因
diffup <-  All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC > foldChange)),]
write.csv(diffup, "diffup.csv")
#
diffdown <- All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig < -foldChange)),]
write.csv(diffdown, "diffdown.csv")

## 导入R包
library(ggplot2)
library(ggrepel)
##  绘制火山图
## 进行分类别
logFC <- resOrdered$log2FoldChange
deg.padj <- resOrdered$pvalue
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > 0.05 | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "Not"
data$group[(data$padj <= 0.05 & data$logFC > 2)] <-  "Up"
data$group[(data$padj <= 0.05 & data$logFC < -2)] <- "Down"
x_lim <- max(logFC,-logFC)

label = subset(resOrdered,pvalue <0.05 & abs(log2FoldChange) > 2)
label1 = rownames(label)

colnames(resOrdered)[1] = 'baseMean'
Significant=ifelse((resOrdered$pvalue < 0.05 & abs(resOrdered$log2FoldChange)> 2), ifelse(resOrdered$log2FoldChange > 2,"Up","Down"), "Not")

ggplot(resOrdered, aes(log2FoldChange, -log10(pvalue)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))

#####GESA method
library(fgsea)          # GSEA分析主程序
library(data.table)     # 数据处理
library(ggplot2)        # 画图处理
library(dplyr)          # 数据处理
library(msigdb)         # 包含基因集合，通常和GSEA分析共同使用
library(GSEABase)       # 可以提供GSEA基础结构和函数,也会被其他包调用
library(stringi)##加载包
library(org.Hs.eg.db)
library(clusterProfiler)
# 查看org.Hs.eg.db 包提供的转换类型
keytypes(org.Hs.eg.db)

# 采用bitr()函数进行转换
library(fgsea)
library(stats)
library(BiocParallel)
library(enrichplot)

rownames(All_diffSig)=stri_sub(rownames(All_diffSig),1,15)##保留前15位
Ensembl_ID <- rownames(All_diffSig)
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
# 查看转换的结果
head(gene_symbol)
All_diffSig.genesymol<-All_diffSig[gene_symbol$ENSEMBL,]
which(gene_symbol$ENTREZID=="641367")#509,719
which(gene_symbol$ENTREZID=="653145")#437,797
gene_symbol<-gene_symbol[-c(719,797),]
All_diffSig.genesymol<-All_diffSig.genesymol[-c(719,797),]
rownames(All_diffSig.genesymol)<-gene_symbol$ENTREZID
ge = All_diffSig.genesymol$log2FoldChange
names(ge) = rownames(All_diffSig.genesymol)
ge = sort(ge,decreasing = T)
BiocParallel::register(BiocParallel::SerialParam())#the first time before running fgsea. This will disable parallel computation.
#GSEA分析——GO
Go_gseresult <- gseGO(ge, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA分析——KEGG
KEGG_gseresult <- gseKEGG(ge,  minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA分析——Reactome
library(ReactomePA)
Go_Reactomeresult <- gsePathway(ge, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
write.table (Go_gseresult, file ="Go_gseresult1.csv", sep =",", row.names =TRUE)
write.table (KEGG_gseresult, file ="KEGG_gseresult.csv", sep =",", row.names =TRUE)
write.table (Go_Reactomeresult, file ="Go_Reactomeresult.csv", sep =",", row.names =TRUE)

#波浪图
ridgeplot(Go_gseresult,10)
ridgeplot(KEGG_gseresult,10)#no
ridgeplot(Go_Reactomeresult,10) #输出前十个结果
#gseaplot2还可以同时显示复数个功能组的富集曲线，并标记P值：
gseaplot2(KEGG_gseresult, 1:6, pvalue_table = TRUE)
###############GSVA
#KEGG
KEGG_df_all <-  msigdbr(species = "Homo sapiens", # Homo sapiens or Mus musculus
                        category = "C2",
                        subcategory = "CP:KEGG") 
KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)

kegg_list = split(KEGG_df$gene_symbol,KEGG_df$gs_exact_source)
lapply(kegg_list[1:3], head)
library(GSVA)
library(msigdbr)
grep("SLC35E2",rownames(mRNA))
mRNA<-mRNA[-16272,]
mRNA<-as.matrix(mRNA)
rownames(mRNA)<-str_replace_all(rownames(mRNA),"\\|\\S+","")
KEGG_ES <- gsva(expr=mRNA, 
                gset.idx.list=kegg_list, 
                parallel.sz=2) #自己电脑parallel.sz写5就好，线程数
#immune
immune_marker<-read.csv("X:/immuneMarker/marker.csv")
immune_ES <- gsva(expr=mRNA, 
                gset.idx.list=immune_marker, 
                parallel.sz=2) #自己电脑parallel.sz写5就好，线程数

#GO

GO_df_all <- msigdbr(species = "Homo sapiens",
                     category = "C5")  
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
go_list <- split(GO_df$gene_symbol, GO_df$gs_name) ##按照gs_name给gene_symbol分组


go_list = split(GO_df$gene_symbol,GO_df$gs_exact_source)
lapply(go_list[1:3], head)
GO_ES <- gsva(expr=mRNA, 
              gset.idx.list=go_list, 
              parallel.sz=32) #自己电脑parallel.sz写5就好，线程数
rownames(annnote332)<-str_replace_all(rownames(annnote332),"\\-","\\.")
GO_ES332<-GO_ES[,rownames(annnote332)]
KEGG_ES332<-KEGG_ES[,rownames(annnote332)]
immune_ES332<-immune_ES[,rownames(annnote332)]
###GO 通路和注释信息
library(GO.db)
goterms <- Term(GOTERM)
GOlist=as.data.frame(goterms)
rownames(GO_ES332)<-GOlist[rownames(GO_ES332),]
###KEGG ID和注释信息
KEGG_term<-unique(KEGG_df$gs_name)
KEGG_ID<-unique(KEGG_df$gs_exact_source)
KEGG.annote<-data.frame(ID=KEGG_ID,term=KEGG_term)
rownames(KEGG.annote)<-KEGG.annote$ID
rownames(KEGG_ES332)<-KEGG.annote[rownames(KEGG_ES332),]$term

#### 进行limma差异处理 ####
##设定 low/high
library(limma)
group_list <-as.factor(annnote332$class) 
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- rownames(annnote332)
contrast.matrix <- makeContrasts(high-low,  #"high/low"
                                 levels = design)

fit1 <- lmFit(KEGG_ES332,design)                 #拟合模型
fit1 <- lmFit(immune_ES332,design)                 #拟合模型

fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正

summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
tempOutput <- topTable(efit,n = Inf, adjust = "fdr")
degs <- na.omit(tempOutput) 
write.csv(degs,"gsva_go_degs.results.csv")
#### 对GSVA的差异分析结果进行热图可视化 #### 
##设置筛选阈值
padj_cutoff=0.05
log2FC_cutoff=log2(1)
##immune
immune.sig<-degs[degs$adj.P.Val < padj_cutoff & abs(degs$logFC)>log2FC_cutoff, ]
immune.sig<-immune.sig[order(abs(immune.sig$logFC),decreasing=T ),]
keepimmune <- rownames(immune.sig)
#KEGG
KEGG.sig<-degs[degs$adj.P.Val < padj_cutoff & abs(degs$logFC)>log2FC_cutoff, ]
KEGG.sig<-KEGG.sig[order(abs(KEGG.sig$logFC),decreasing=T ),]
keepKEGG <- rownames(KEGG.sig)
#GO
GO.sig<-degs[degs$adj.P.Val < padj_cutoff & abs(degs$logFC)>log2FC_cutoff, ]
GO.sig<-GO.sig[order(abs(GO.sig$logFC),decreasing=T),]
keepGO <- GO.sig$ID

length(keepimmune)
dat <- immune_ES332[keepimmune,] #选取前30进行展示KEGG_ES332[keepKEGG[1:30],]
win.graph(width=4.875, height=2.5,pointsize=8)
annnote332class<-data.frame(class=annnote332$class)
rownames(annnote332class)<-colnames(dat)
#rownames(dat)[3]<-"oxidoreductase activity, acting on paired donors"
pheatmap::pheatmap(dat, 
                   fontsize_row = 8,
                   height = 10,
                   width=18,
                   cluster_cols = FALSE,
                   annotation_col = annnote332class,
                   show_colnames = F,
                   show_rownames = T
                )
#####immu infiltration analysis
#timer
mRNA332<-mRNA[,rownames(annnote332)]
table(apply(mRNA332, 1, function(x){sum(x==0)<100}))
tmp <- apply(mRNA332, 1, function(x){sum(x==0)<100})
mRNA332 <- mRNA332[tmp,]
write.csv(mRNA332,file = "./mRNA332.csv",sep = "\t")
Timer_score<-read.csv("./Timer_score_matrix.csv")
#Timer_score$class<-annnote332[Timer_score$sampleID,]$class
Timer_score<-t(Timer_score)
Timer_score<-as.data.frame(Timer_score)
colnames(Timer_score)<-Timer_score[1,]
Timer_score<-Timer_score[-1,]
Timer_score$Immunecell<-rownames(Timer_score)
boxplot.df.timer<-melt(Timer_score,id.vars='Immunecell')
annnote332$variable<-rownames(annnote332)
boxplot.df.timer.class<-merge(boxplot.df.timer,annnote332,by=c("variable"="variable"))
boxplot.df.timer.class<-boxplot.df.timer.class[,-4]
boxplot.df.timer.class$value<-as.numeric(boxplot.df.timer.class$value)
df_p_val <- boxplot.df.timer.class %>% 
  group_by(Immunecell) %>% 
  wilcox_test(formula = value ~ class) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='immunecell')
#p1 <- ggplot(boxplot.df,aes(x=miName,y=logvalue,fill=Class))+
#  geom_boxplot(width=0.6,alpha=0.8)
#win.graph(width=4.875, height=2.5,pointsize=8)
ggplot()+
  geom_boxplot(data = boxplot.df.timer.class,mapping = aes(x=Immunecell,y=value,fill=class),width=0.5)+
  scale_fill_manual(values = c('#E21C21','#3A7CB5'))+
  stat_pvalue_manual(
    df_p_val, x = "Immunecell", y.position = df_p_val$y.position+0.2,
    label = "p.signif",label.size=3,
    position = position_dodge(0.8)
  )+
  labs(y='Value',x='')+
  theme_test()+
  theme(axis.text = element_text(color = 'black',size = 15),#????????????С
        plot.caption = element_markdown(face = 'bold'),
        legend.position = c(0.85,0.95),legend.text = element_text (size = 15),#ͼ????????С
        legend.title = element_text(size=15),#ͼ???????ִ?С
        axis.title.y=element_text(size=15),#y??label??С
        legend.direction = 'horizontal')+#ͼ??????
  theme (axis.text.x = element_text (angle=45, hjust=1, vjust=1))
####cibersort
mRNA332<-as.matrix(mRNA332)
#mRNA322.new<-data.frame(Symbol=rownames(mRNA332),mRNA332)
library('CIBERSORT')
#读取包自带的LM22文件（免疫细胞特征基因文件）
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
results.cb <- cibersort(sig_matrix, mRNA332,perm = 1000,QN = T)
results.cb<-as.data.frame(results.cb)
results<-results.cb[,1:22]
results<-t(results)
results<-as.data.frame(results)
results$Cells<-rownames(results)
boxplot.df.ciber<-melt(results,id.vars='Cells')
boxplot.df.ciber.class<-merge(boxplot.df.ciber,annnote332,by=c("variable"="variable"))

boxplot.df.ciber.class<-boxplot.df.ciber.class[,-4]
boxplot.df.ciber.class$value<-as.numeric(boxplot.df.ciber.class$value)
df_p_val <- boxplot.df.ciber.class %>% 
  group_by(Cells) %>% 
  wilcox_test(formula = value ~ class) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Cells')
#p1 <- ggplot(boxplot.df,aes(x=miName,y=logvalue,fill=Class))+
#  geom_boxplot(width=0.6,alpha=0.8)
#win.graph(width=4.875, height=2.5,pointsize=8)
ggplot()+
  geom_boxplot(data = boxplot.df.ciber.class,mapping = aes(x=Cells,y=value,fill=class),width=0.5)+
  scale_fill_manual(values = c('#E21C21','#3A7CB5'))+
  stat_pvalue_manual(
    df_p_val, x = "Cells", y.position = df_p_val$y.position+0.2,
    label = "p.signif",label.size=3,
    position = position_dodge(0.8)
  )+
  labs(y='Value',x='')+
  theme_test()+
  theme(axis.text = element_text(color = 'black',size = 15),#????????????С
        plot.caption = element_markdown(face = 'bold'),
        legend.position = c(0.85,0.95),legend.text = element_text (size = 15),#ͼ????????С
        legend.title = element_text(size=15),#ͼ???????ִ?С
        axis.title.y=element_text(size=15),#y??label??С
        legend.direction = 'horizontal')+#ͼ??????
  theme (axis.text.x = element_text (angle=45, hjust=1, vjust=1))



