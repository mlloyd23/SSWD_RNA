library(dplyr)
setwd("C:/Users/Melanie/Documents/R/deseq/collagen enrichment")

#phenotype
pheno<-read.table("pheno_results_no5.txt", sep="\t", header=TRUE)
head(pheno)
pheno<-pheno[complete.cases(pheno), ]
dim(pheno)

pheno$abs.stat<-abs(pheno$stat)
hist(pheno$abs.stat)

t.test(abs.stat~collagen, data=pheno)
wilcox.test(abs.stat~collagen, data=pheno)
boxplot(abs.stat~collagen, data=pheno)

boxplot(stat~collagen, data=pheno)
hist(pheno$stat)
qqnorm(pheno$stat)
qqline(pheno$stat)

#density plots
collagen<-subset(pheno, collagen=="1")
dim(collagen)
not<-subset(pheno, collagen=="0")
dim(not)
  
plot(density(collagen$stat), ylim=c(0,0.6))
lines(density(not$stat))


library(ggplot2)
pheno$collagen<-as.factor(pheno$collagen)
ggplot(pheno, aes(x=stat, fill=collagen)) +
  geom_density(alpha=0.4)+
  scale_fill_manual(values=c("#999999", "red3"), 
                    labels=c("Non-collagen\ngenes", "Collagen\nrelated genes"))+
  xlim(-5,5)+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))



##immune enrichment???
immune.pheno<-read.table("immune_enrich_pheno_no5.txt", header=TRUE)
head(immune.pheno)
immune.pheno<-immune.pheno[complete.cases(immune.pheno), ]
dim(immune.pheno)
immune<-subset(immune.pheno, immune=="1")
not.immune<-subset(immune.pheno, immune=="0")
plot(density(immune$stat), ylim=c(0,0.6))
lines(density(not.immune$stat))

t.test(stat~immune, data=immune.pheno)

#anova interaction
aov<-read.table("anova_int.txt", header=TRUE)
aov<-aov[complete.cases(aov), ]
aov$abs.stat<-abs(aov$stat)

hist(aov$abs.stat)

t.test(abs.stat~collagen, data=aov)
wilcox.test(abs.stat~collagen,data=aov)

#density plots
collagen.aov<-subset(aov, collagen=="1")
dim(collagen.aov)
not.aov<-subset(aov, collagen=="0")
dim(not.aov)

plot(density(collagen.aov$stat), ylim=c(0,2))
lines(density(not.aov$stat))


#pn10
pn10<-read.table("pn10.txt", header=TRUE)
dim(pn10)
pn10<-pn10[complete.cases(pn10), ]
hist(pn10$stat)
pn10$abs.stat<-abs(pn10$stat)
head(pn10)

t.test(abs.stat~collagen, data=pn10)
wilcox.test(abs.stat~collagen, data=pn10)
boxplot(abs.stat~collagen, data=pn10)

#density plots
collagen.pn10<-subset(pn10, collagen=="1")
dim(collagen.pn10)
not.pn10<-subset(pn10, collagen=="0")
dim(not.pn10)

plot(density(collagen.pn10$stat), ylim=c(0,0.6))
lines(density(not.pn10$stat))

#pn21
pn21<-read.table("pn21.txt", header=TRUE)
hist(pn21$stat)
pn21$abs.stat<-abs(pn21$stat)
head(pn21)

t.test(abs.stat~collagen, data=pn21)
wilcox.test(abs.stat~collagen, data=pn21)
boxplot(abs.stat~collagen, data=pn21)


#pn20
pn20<-read.table("pn20.txt", header=TRUE)
hist(pn20$stat)
pn20$abs.stat<-abs(pn20$stat)

t.test(abs.stat~collagen, data=pn20)
wilcox.test(abs.stat~collagen, data=pn20)
boxplot(abs.stat~collagen,data=pn20)
