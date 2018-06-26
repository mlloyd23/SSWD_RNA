library(CHNOSZ)
library(DESeq2)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(tidyr)
library(RColorBrewer)


blast_pm<-read.blast("tblastx_pm.txt", evalue = 1e-5, max.hits = 1, 
                     min.length = NA, quiet = FALSE)

data<-read.table("grouped-counts.txt", header=TRUE)
colnames(data)<-c("03-5-08", "03-5-11", "07-5-08", "07-5-11", "08-5-08", 
                  "08-5-11", "08-5-14", "08-5-17", "08-5-20", "09-5-08", 
                  "09-5-14", "09-5-17", "09-5-20", "10-5-08", "10-5-11", 
                  "10-5-14", "10-5-17", "10-5-20", "14-5-08", "14-5-11", 
                  "15-5-08", "15-5-11", "15-5-14", "15-5-17", "15-5-20", 
                  "19-5-11", "19-5-14", "19-5-17", "19-5-20", "20-5-08", 
                  "20-5-11", "20-5-14", "20-5-17", "20-5-20", "22-5-08", 
                  "22-5-11", "22-5-14", "23-5-17", "23-5-20", "24-5-08", 
                  "24-5-11", "24-5-14", "24-5-17", "24-5-20", "26-5-08", 
                  "26-5-11", "27-5-08", "27-5-11", "27-5-14", "27-5-17", 
                  "27-5-20", "28-5-08", "28-5-11", "28-5-14", "28-5-17", 
                  "29-5-08", "29-5-11", "29-5-14", "31-6-12", "31-6-15", 
                  "31-6-18", "31-6-21", "31-6-24", "32-6-12", "32-6-15", 
                  "32-6-18", "32-6-21", "33-6-12", "33-6-15", "33-6-18", 
                  "33-6-21", "33-6-24", "34-6-12", "34-6-15", "34-6-18", 
                  "34-6-21", "34-6-24", "35-6-12", "35-6-15", "35-6-18", 
                  "35-6-21", "36-6-12", "36-6-15", "36-6-18", "37-6-12", 
                  "37-6-15", "37-6-18", "37-6-21", "38-6-12", "38-6-15", 
                  "38-6-18", "38-6-21", "38-6-24")

map<-read.table("map.txt", header=TRUE)
rownames(map)<-c("03-5-08", "03-5-11", "07-5-08", "07-5-11", "08-5-08", 
                 "08-5-11", "08-5-14", "08-5-17", "08-5-20", "09-5-08", 
                 "09-5-14", "09-5-17", "09-5-20", "10-5-08", "10-5-11", 
                 "10-5-14", "10-5-17", "10-5-20", "14-5-08", "14-5-11", 
                 "15-5-08", "15-5-11", "15-5-14", "15-5-17", "15-5-20", 
                 "19-5-11", "19-5-14", "19-5-17", "19-5-20", "20-5-08", 
                 "20-5-11", "20-5-14", "20-5-17", "20-5-20", "22-5-08", 
                 "22-5-11", "22-5-14", "23-5-17", "23-5-20", "24-5-08", 
                 "24-5-11", "24-5-14", "24-5-17", "24-5-20", "26-5-08", 
                 "26-5-11", "27-5-08", "27-5-11", "27-5-14", "27-5-17", 
                 "27-5-20", "28-5-08", "28-5-11", "28-5-14", "28-5-17", 
                 "29-5-08", "29-5-11", "29-5-14", "31-6-12", "31-6-15", 
                 "31-6-18", "31-6-21", "31-6-24", "32-6-12", "32-6-15", 
                 "32-6-18", "32-6-21", "33-6-12", "33-6-15", "33-6-18", 
                 "33-6-21", "33-6-24", "34-6-12", "34-6-15", "34-6-18", 
                 "34-6-21", "34-6-24", "35-6-12", "35-6-15", "35-6-18", 
                 "35-6-21", "36-6-12", "36-6-15", "36-6-18", "37-6-12", 
                 "37-6-15", "37-6-18", "37-6-21", "38-6-12", "38-6-15", 
                 "38-6-18", "38-6-21", "38-6-24")

map$Day<-as.factor(map$Day)
map$individual<-as.factor(map$individual)

#Pmhits
Pm<-blast_pm[,1]
df_Pm<-cbind(Row.Names=rownames(data),data)
df_Pm<-subset(df_Pm, df_Pm$Row.Names %in% Pm)
rownames(df_Pm)<-df_Pm$Row.Names
df_Pm <- df_Pm[ -c(1) ]

##Phenotype
countsTable <- DESeqDataSetFromMatrix(
  countData = df_Pm,
  colData = map,
  design = ~  individual + Phenotype)

de<-DESeq(countsTable)
pheno_results<-results(de, contrast=c("Phenotype","Sick","Healthy"))
head(pheno_results)
summary(pheno_results)
sigtab <- pheno_results[which(pheno_results$padj < 0.1), ]

write.table(pheno_results, "pheno_results.txt", sep="\t")

##Phenotype, no dead samples
countsTable_sub<-countsTable[,colData(countsTable)$Pheno_num!=5]
de_no5<-DESeq(countsTable_sub)
pheno_results_no5<-results(de_no5, contrast=c("Phenotype","Sick","Healthy"))
head(pheno_results_no5)
summary(pheno_results_no5)
write.table(pheno_results_no5, "pheno_results_no5_Rupdated.txt", sep="\t")

###Phenotype number
map$pn<-as.factor(map$pn)
countsTable_pn <- DESeqDataSetFromMatrix(
  countData = df_Pm,
  colData = map,
  design = ~  individual + pn)

dds_pn <- DESeq(countsTable_pn) 
res_pn21 <- results(dds_pn, contrast=c("pn","2","1"))
head(res_pn21)
summary(res_pn21)
sigtab <- res_pn21[which(res_pn21$padj < 0.1), ]
sigtab

res_pn10 <- results(dds_pn, contrast=c("pn","1","0"))
head(res_pn10)
summary(res_pn10)
sigtab <- res_pn10[which(res_pn10$padj < 0.1), ]
sigtab

res_pn20 <- results(dds_pn, contrast=c("pn","2","0"))
head(res_pn20)
summary(res_pn20)

write.table(res_pn20, "pn20.txt", sep="\t")

##ANOVA

anova_map<-read.table("ANOVA_map.txt", header=TRUE)
anova_map$Day<-as.factor(anova_map$Day)
anova_map$Individual<-as.factor(anova_map$Individual)

anova_samples<-as.character(anova_map[,1])
anova_df<- select(df_Pm, one_of(anova_samples))

countsTable_anova <- DESeqDataSetFromMatrix(
  countData = anova_df,
  colData = anova_map,
  design = ~  Time + Final_phenotype + Time:Final_phenotype)

int_anova <- DESeq(countsTable_anova, test="LRT", reduced = ~ Time + Final_phenotype) 
res_int <- results(int_anova)
head(res_int)
summary(res_int)
sigtab <- res_int[which(res_int$padj < 0.1), ]
write.table(res_int,"ANOVA_int.txt",sep="\t")

##FIGURE 1 BOXPLOTS
#Fig.1A
goi <- c("Cluster-33358.68302", "Cluster-33358.54259",
         "Cluster-33358.50323")

tcounts <- t(log2((counts(dds_pn[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds_pn), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

tcounts$Day<-as.factor(tcounts$Day)
tcounts$pn<-as.factor(tcounts$pn)

tcounts %>% 
  select(Row.names, Day, pn,Phenotype, gene, expression) %>% 
  head %>% 
  knitr::kable()

palette<-c("darkorange", "darkorchid4")

labels <- c("Cluster-33358.68302"= "complement \ncomponent C3", 
            "Cluster-33358.54259"= "lysozyme",
            "Cluster-33358.50323"="ras-related \nprotein Rab-18-B")

title<-"Immune genes"
ggplot(tcounts, aes(Phenotype, expression, fill=Phenotype)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y",  labeller=labeller(gene = labels)) + 
  labs(x=NULL, 
       y="Expression (log normalized counts)", 
       fill="Symptom Stage")+
  theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title = element_text(colour="#666666", size=12, face="bold"))+
  scale_fill_manual(values= palette, labels=c("Healthy", "Sick"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  ggtitle(title)


#Fig.1B
goi <- c("Cluster-33358.64315", "Cluster-33358.57839",
         "Cluster-33358.54397")

tcounts <- t(log2((counts(dds_pn[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds_pn), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

tcounts$Day<-as.factor(tcounts$Day)
tcounts$pn<-as.factor(tcounts$pn)

tcounts %>% 
  select(Row.names, Day, pn,Phenotype, gene, expression) %>% 
  head %>% 
  knitr::kable()

palette<-c("darkorange", "darkorchid4")

labels <- c("Cluster-33358.64315"= " translation\ninitiation\nfactor 3", 
            "Cluster-33358.57839"= "RNA-binding\ncabeza isoform",
            "Cluster-33358.54397"="activated RNA\npolymerase II\ntranscriptional\ncoactivator")

title<-"RNA genes"
ggplot(tcounts, aes(Phenotype, expression, fill=Phenotype)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y",  labeller=labeller(gene = labels)) + 
  labs(x=NULL, 
       y="Expression (log normalized counts)", 
       fill="Symptom Stage")+
  theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title = element_text(colour="#666666", size=12, face="bold"))+
  scale_fill_manual(values= palette, labels=c("Healthy", "Sick"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  ggtitle(title)



##FIGURE 2a PCA PLOT
 
rld <- varianceStabilizingTransformation(countsTable_pn, blind=FALSE)


p=plotPCA(rld, intgroup=c("pn"))
p1= p + geom_point(size = 2, alpha = 0.7) + 
  scale_colour_manual(values = c("#fecc5c", "#7fcdbb", "#2c7fb8"),
                      name="Disease stage",
                      labels=c("Healthy","Early stage","Late stage"))+
  theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title = element_text(colour="#666666", size=12, face="bold"))+
  scale_fill_manual(values= palette, labels=c("Healthy", "Sick"))+
  ylim(-33,33)+
  xlim(-33,33)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p1+stat_ellipse(type = "t")

p1