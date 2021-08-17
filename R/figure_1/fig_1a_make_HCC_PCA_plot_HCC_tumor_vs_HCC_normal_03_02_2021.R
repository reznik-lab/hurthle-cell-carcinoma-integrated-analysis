library(ggplot2)
library(ggrepel)

####
# load protected additional clinical information
####

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/PHI"
fileName<-"clin.Rda"
fileName<-file.path(filePath,fileName)
load(fileName)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/PHI"
fileName<-"clingen.Rda"
fileName<-file.path(filePath,fileName)
load(fileName)

######

# add recurrence to clin
clin$Histology2 = clin$HISTOLOGY
recur = rownames(clingen)[which(clingen$Progression == "Y")]
clin[recur,'Histology2'] = 'HWIDE Recurrent'

clin_HWIDE<-rownames(clin[clin$Histology2 %in% "HWIDE",])
clin_HWIDE_recurrent<-rownames(clin[clin$Histology2 %in% "HWIDE Recurrent",])

# ADD YD53 and YD62 after EMR review
clin_HWIDE_recurrent<-c("YD2","YD20","YD22","YD31","YD39","YD94","YD53","YD62")


######

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/metabotools_HCC"
#fileName<-"HCC_met_prefiltered_PQN_normalized.Rda"
fileName<-"HCC_met_Elisa_PNQ_normalized_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

# met 737 metabolite x 75 samples
# dim(met)
load(fileName)

#####
# load metabolite info
#####

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_metabolite_info_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

load(file=fileName)

#####
# load smaple info
#####
filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_sample_info_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

# sampinfo data frame
load(fileName)

sampleInfo<-t(sampinfo)
sampleInfo<-as.data.frame(sampleInfo,stringsAsFactors = FALSE)

# fix label
sampleInfo[sampleInfo$CLIENT_IDENTIFIER %in% "YD53",]$HISTOLOGY<-"widely invasive HCC"
sampleInfo[sampleInfo$CLIENT_IDENTIFIER %in% "YD35",]$HISTOLOGY<-"normal"
sampleInfo[sampleInfo$HISTOLOGY %in% "Normal",]$HISTOLOGY<-"normal"

# add histology 2 for ggplot labeling
sampleInfo$HISTOLOGY2<-sampleInfo$HISTOLOGY
sampleInfo[sampleInfo$HISTOLOGY2 %in% "widely invasive HCC",]$HISTOLOGY2<-"HWIDE"
sampleInfo[sampleInfo$CLIENT_IDENTIFIER %in% clin_HWIDE_recurrent,]$HISTOLOGY2<-"HWIDE Recurrent"
sampleInfo[sampleInfo$HISTOLOGY2 %in% "minimally invasive HCC",]$HISTOLOGY2<-"HMIN"
sampleInfo[sampleInfo$HISTOLOGY2 %in% "HA",]$HISTOLOGY2<-"HA"
sampleInfo[sampleInfo$HISTOLOGY2 %in% "Poorly diff",]$HISTOLOGY2<-"PD"
sampleInfo[sampleInfo$HISTOLOGY2 %in% "PTCTV",]$HISTOLOGY2<-"PTCTV"
sampleInfo[sampleInfo$HISTOLOGY2 %in% "normal",]$HISTOLOGY2<-"Normal"



# from "TUMOR" TISSUE
# HA, 5
# widely invasive HCC (HWIDE/HCC), 18 
# HCC minimally invasive (HMIN/HCC), 14
# HCC, 1 (HWIDE)           
# Poorly diff (PD), 8
# PTCTV, 4 

# from "NORMAL" TISSUE
# normal, 27


tumorSampleName<-sampleInfo[sampleInfo$TISSUE %in% c("TUMOR"),]$CLIENT_IDENTIFIER
normalSampleName<-sampleInfo[sampleInfo$TISSUE %in% c("NORMAL"),]$CLIENT_IDENTIFIER


sampleInfo_Tumor<-sampleInfo[sampleInfo$TISSUE %in% c("TUMOR"),]
sampleInfo_Normal<-sampleInfo[sampleInfo$TISSUE %in% c("NORMAL"),]


hcc_group<-c("HCC","minimally invasive HCC", "widely invasive HCC")
pd_group<-c("Poorly diff")
ptctv_group<-c("PTCTV")
ha_group<-c("HA")

hcc_HWIDE_group<-c("widely invasive HCC")
hcc_HMIN_group<-c("minimally invasive HCC")
pd_and_ptctv_group<-c("Poorly diff","PTCTV")

# drop YD29 (HCC HMIN), YD65 (PDTC) tumor samples at the beginning
# drop YD30, YD66 matched normal samples for analysis    
removedSamples<-c("YD29","YD65")

tumorSampleName_hcc<-sampleInfo_Tumor[sampleInfo_Tumor$HISTOLOGY %in% hcc_group,]$CLIENT_IDENTIFIER
tumorSampleName_hcc<-tumorSampleName_hcc[!( tumorSampleName_hcc %in% removedSamples)]

tumorSampleName_pd<-sampleInfo_Tumor[sampleInfo_Tumor$HISTOLOGY %in% pd_group,]$CLIENT_IDENTIFIER
tumorSampleName_pd<-tumorSampleName_pd[!( tumorSampleName_pd %in% removedSamples)]


tumorSampleName_ptctv<-sampleInfo_Tumor[sampleInfo_Tumor$HISTOLOGY %in% ptctv_group,]$CLIENT_IDENTIFIER
tumorSampleName_ha<-sampleInfo_Tumor[sampleInfo_Tumor$HISTOLOGY %in% ha_group,]$CLIENT_IDENTIFIER


tumorSampleName_hcc_HWIDE<-sampleInfo_Tumor[sampleInfo_Tumor$HISTOLOGY %in% hcc_HWIDE_group,]$CLIENT_IDENTIFIER

tumorSampleName_hcc_HWIDE_no_recurrent<-tumorSampleName_hcc_HWIDE[!(tumorSampleName_hcc_HWIDE %in% clin_HWIDE_recurrent)]
tumorSampleName_hcc_HWIDE_recurrent<-tumorSampleName_hcc_HWIDE[(tumorSampleName_hcc_HWIDE %in% clin_HWIDE_recurrent)]

tumorSampleName_hcc_HMIN<-sampleInfo_Tumor[sampleInfo_Tumor$HISTOLOGY %in% hcc_HMIN_group,]$CLIENT_IDENTIFIER
tumorSampleName_hcc_HMIN<-tumorSampleName_hcc_HMIN[!( tumorSampleName_hcc_HMIN %in% removedSamples)]

tumorSampleName_pd_and_ptctv<-sampleInfo_Tumor[sampleInfo_Tumor$HISTOLOGY %in% pd_and_ptctv_group,]$CLIENT_IDENTIFIER
tumorSampleName_pd_and_ptctv<-tumorSampleName_pd_and_ptctv[!( tumorSampleName_pd_and_ptctv %in% removedSamples)]

# 16 samples in HCC with matched tumor and normal samples
hcc_normal_group<-c("NT_HW_Mut","NT_HW_WT","NT_HMin_Mut","NT_HMin_WT")
normalSampleName_hcc<-sampleInfo_Normal[sampleInfo_Normal$Group %in% hcc_normal_group,]$CLIENT_IDENTIFIER

matchedNormalSampleName_hcc<-sampleInfo_Normal[sampleInfo_Normal$Group %in% hcc_normal_group,]$CLIENT_IDENTIFIER
matchedTumorSampleName_hcc<-sampleInfo_Normal[sampleInfo_Normal$Group %in% hcc_normal_group,]$MATCHED

# drop YD29 (HCC HMIN), YD65 (PDTC) tumor samples at the beginning
# drop YD30, YD66 matched normal samples for analysis    
removedSamples<-c("YD29","YD65","YD30","YD66")
matchedNormalSampleName_hcc<-matchedNormalSampleName_hcc[!( matchedNormalSampleName_hcc %in% removedSamples)]
matchedTumorSampleName_hcc<-matchedTumorSampleName_hcc[!( matchedTumorSampleName_hcc %in% removedSamples)]

# decide to keep YD30, although YD29 sample has dropped
#normalSampleName_hcc<-normalSampleName_hcc[!( normalSampleName_hcc %in% removedSamples)]


######
### PCA checking
######

selectedSamples<-c(tumorSampleName_hcc,normalSampleName_hcc)

newMet<-log2(met[,selectedSamples])

rowSD<-apply(newMet,1,sd)

#sdCutoff<-quantile(sort(rowSD,decreasing = TRUE),percentileThreshold)
sdCutoff<-0
metaboliteSubset<-names(rowSD[rowSD>sdCutoff])

numberOfmetabolitesLeft<-length(metaboliteSubset)
numberOfmetabolitesTotal<-length(rownames(newMet))

message(sprintf("Use %s / %s metabolites in the data",numberOfmetabolitesLeft,numberOfmetabolitesTotal))

newMetPCA<-newMet[metaboliteSubset,]
dim(newMetPCA)



if(FALSE){
  rowSD<-apply(newMet,1,sd)
  percentileThreshold<-0.8
  sdCutoff<-quantile(sort(rowSD,decreasing = TRUE),percentileThreshold)
  metaboliteSubset<-names(rowSD[rowSD>=sdCutoff])
  
  newMetPCA<-newMet[metaboliteSubset,]
  dim(newMetPCA)
}


if(FALSE){

completePercentage<-apply(newMet, 1, function(x){
  (length(x) - sum(x == min(x)) +1 )/length(x)
})

plotDat<-data.frame(metabolite=names(completePercentage),value=1-completePercentage,stringsAsFactors = FALSE)
plotDat<-plotDat[plotDat$value!=0,]
plotDat<-plotDat[order(-plotDat$value),]
ggplot(plotDat,aes(x=metabolite,y=value,label=metabolite))+geom_point()+geom_text_repel()+theme_classic()

threshold<-0.9
completeMetabolites<-names(completePercentage[completePercentage>threshold])

numberOfmetabolitesLeft<-length(completeMetabolites)
numberOfmetabolitesTotal<-length(rownames(newMet))

message(sprintf("Use %s / %s metabolites in the data",numberOfmetabolitesLeft,numberOfmetabolitesTotal))

metaboliteSet<-completeMetabolites

newMetPCA<-newMet[rownames(newMet) %in% metaboliteSet,]

#newMetPCA<-newMet[rownames(newMet) %in% completeMetabolites,]

message(sprintf("Use %s / %s metabolites for %s samples in the data",nrow(newMetPCA), numberOfmetabolitesTotal, ncol(newMetPCA)))

}


######
# Start of PCA calculation section
# transpose matrix to do sample PCA
######

res.pca = prcomp(t(newMetPCA),scale = TRUE, center = TRUE)

#######
# PCA screeplot
#######

#res.pca<-samppca

#https://stackoverflow.com/questions/60957020/how-to-set-scree-plot-scale-as-same-as-principal-components
#variance explained
varExp = (100*res.pca$sdev^2)/sum(res.pca$sdev^2)

varPC1<-round(varExp[1],digit=2)
varPC2<-round(varExp[2],digit=2)

varDF = data.frame(Dimensions=1:length(varExp),
                   varExp=varExp)

graph <- ggplot(varDF,aes(x=Dimensions,y=varExp)) 
graph <- graph + geom_point() 
graph <- graph + geom_col(fill="steelblue") 
graph <- graph + geom_line() 
graph <- graph + scale_x_continuous(breaks=1:nrow(varDF))
graph <- graph + ylim(c(0,100)) 
graph <- graph + ylab("Variance explained (%)")
graph <- graph + theme_bw()

print(graph)

#####

if(FALSE){
###select PCs that explain at least 80% of the variance

varExp = (res.pca$sdev^2)/sum(res.pca$sdev^2)

cum_var <- cumsum(varExp)
select_pcs <- which(cum_var>0.8)[1]

pre_rank_matrix <- as.matrix(rowSums(abs(res.pca$rotation[,1:select_pcs])))
pre_rank_matrix <- pre_rank_matrix[order(pre_rank_matrix[,1],decreasing = TRUE),]

pre_rank_dat<-data.frame(names(pre_rank_matrix),pre_rank_matrix,stringsAsFactors = FALSE)
colnames(pre_rank_dat)<-c("name","rankValue")

#######
# handle differential metabolite abundance list
#######

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_metabolite_info_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

load(file=fileName)

metinfoSubset<-metinfo[,c("name","H_KEGG")]

pre_rank_dat<-merge(pre_rank_dat,metinfoSubset,by="name",all.x = TRUE)
pre_rank_dat<-pre_rank_dat[order(pre_rank_dat$rankValue,decreasing = TRUE),]
pre_rank_dat<-pre_rank_dat[!(is.na(pre_rank_dat$H_KEGG)),]

# C01885?
pre_rank_dat<-pre_rank_dat[!duplicated(pre_rank_dat$H_KEGG),]


gseaRanks<-pre_rank_dat$rankValue
names(gseaRanks)<-pre_rank_dat$H_KEGG


library(fgsea)

filePath<-"~/work/Ed_lab/HCC_project/data/KEGG_metabolic_pathway/parsed_data"
fileName<-paste("KEGG_metabolic_pathway_list.Rd",sep="")
fileName<-file.path(filePath,fileName)

load(file=fileName)

set.seed(42)


fgseaRes <- fgsea(pathways = metabolitePathway, 
                  stats    = gseaRanks,
                  scoreType = "pos",
                  #eps      = 0.0,
                  minSize  = 5,
                  maxSize  = 200)

fgseaRes<-fgseaRes[order(padj, -abs(NES)), ]

head(fgseaRes[order(padj, -abs(NES)), ], n=10)


fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#####

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG metabolic pathways NES from GSEA") + 
  theme_classic()

}

#####

if(FALSE){
  library(factoextra)
  
  
  fviz_eig(res.pca, addlabels = TRUE)
  fviz_contrib(res.pca,choice="var",axes=1,top=10)
  fviz_contrib(res.pca,choice="var",axes=2,top=10)
  
  fviz_contrib(res.pca, choice = "ind", axes = 1:2)

  fviz_pca_ind(res.pca,repel=TRUE)
  
  
  
  filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/PCA_plot/01_22_2021"
  dir.create(filePath, showWarnings = TRUE, recursive = TRUE)
  
  fileName<-"HCC_PCA_color_by_group_HCC_tumor_vs_HCC_normal_scree_plot_01_22_2021.pdf"
  fileName<-file.path(filePath,fileName)
  
  pdf(fileName,height=5,width=7)
   fviz_eig(res.pca, addlabels = TRUE)
  dev.off()
  
  fileName<-"HCC_PCA_color_by_group_HCC_tumor_vs_HCC_normal_pc1_contribution_01_22_2021.pdf"
  fileName<-file.path(filePath,fileName)
  
  pdf(fileName,height=5,width=7)
    fviz_contrib(res.pca,choice="var",axes=1,top=10)
  dev.off()
  
  
  fileName<-"HCC_PCA_color_by_group_HCC_tumor_vs_HCC_normal_pc2_contribution_01_22_2021.pdf"
  fileName<-file.path(filePath,fileName)
  
  pdf(fileName,height=5,width=7)
    fviz_contrib(res.pca,choice="var",axes=2,top=10)
  dev.off()
  
  #fviz_contrib(res.pca, choice = "ind", axes = 1:2)
  
  
  #fviz_pca_ind(res.pca,repel=TRUE)

  
}

####

#mixed_sampleInfo<-tubeSampleInfo
#rownames(mixed_sampleInfo)<-mixed_sampleInfo$tubeSample

#samppca = prcomp(na.omit(log2(met)),scale. = TRUE, center = TRUE)
samppca = data.frame(res.pca$x[,1:2])
#samppca$Tissue = clin[rownames(samppca),'FORMAT']
#samppca$histology = clin[rownames(samppca),'HISTOLOGY']
#samppca$group = mixed_sampleInfo[rownames(samppca),'group']
samppca$tissue = sampleInfo[rownames(samppca),'TISSUE']
samppca$histology = sampleInfo[rownames(samppca),'HISTOLOGY']
samppca$group = sampleInfo[rownames(samppca),'HISTOLOGY2']

sampleLabel<-sampleInfo[rownames(samppca),]$CLIENT_IDENTIFIER
#batchLabel<-mixed_sampleInfo[rownames(samppca),]$LCMSbatchDate


#ggplot(samppca,aes(PC1,PC2,shape = Tissue)) + geom_point() + theme_minimal(base_size = 16)

theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))

xRange<-ceiling(max(abs(samppca$PC1)))
yRange<-ceiling(max(abs(samppca$PC2)))

#p <- ggplot(samppca,aes(x=PC1, y=PC2, color=Tissue))
#p <- ggplot(samppca,aes(x=PC1, y=PC2, color=group,label=rownames(samppca)))
p <- ggplot(samppca,aes(x=PC1, y=PC2, color=group,label=sampleLabel))
#p <- p + geom_point(aes(shape=mutation))
p <- p + geom_point(aes(shape=tissue))
p <- p + xlab(paste("PC1 (",varPC1,"%)",sep=""))
p <- p + ylab(paste("PC2 (",varPC2,"%)",sep=""))
p <- p + xlim(-xRange,xRange)
p <- p + ylim(-yRange,yRange)
p <- p + geom_hline(yintercept = 0, linetype="dashed",color="grey")
p <- p + geom_vline(xintercept = 0, linetype="dashed",color="grey")
p <- p + geom_text_repel(size=3)
p <- p + theme

print(p)


#####

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/PCA_plot/01_22_2021"
dir.create(filePath, showWarnings = TRUE, recursive = TRUE)

fileName<-"HCC_PCA_color_by_group_HCC_tumor_vs_HCC_normal_01_22_2021.pdf"

fileName<-file.path(filePath,fileName)

pdf(fileName,height=5,width=7)
print(p)
dev.off()

#####

samppca = data.frame(res.pca$x[,1:2])
#samppca$Tissue = clin[rownames(samppca),'FORMAT']
#samppca$histology = clin[rownames(samppca),'HISTOLOGY']
#samppca$group = mixed_sampleInfo[rownames(samppca),'group']
samppca$tissue = sampleInfo[rownames(samppca),'TISSUE']
samppca$histology = sampleInfo[rownames(samppca),'HISTOLOGY']
samppca$group = sampleInfo[rownames(samppca),'HISTOLOGY2']

sampleLabel<-sampleInfo[rownames(samppca),]$CLIENT_IDENTIFIER
#batchLabel<-mixed_sampleInfo[rownames(samppca),]$LCMSbatchDate


#ggplot(samppca,aes(PC1,PC2,shape = Tissue)) + geom_point() + theme_minimal(base_size = 16)

theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))

xRange<-ceiling(max(abs(samppca$PC1)))
yRange<-ceiling(max(abs(samppca$PC2)))

#p <- ggplot(samppca,aes(x=PC1, y=PC2, color=Tissue))
#p <- ggplot(samppca,aes(x=PC1, y=PC2, color=group,label=rownames(samppca)))
p <- ggplot(samppca,aes(x=PC1, y=PC2, color=group))
p <- p + scale_color_manual(values=c("steelblue1","tomato","mediumpurple1","grey55"))
#p <- p + geom_point(aes(shape=mutation))
p <- p + geom_hline(yintercept = 0, linetype="dashed",color="grey85")
p <- p + geom_vline(xintercept = 0, linetype="dashed",color="grey85")
p <- p + geom_point(aes(shape=tissue),size=1.5)
p <- p + xlab(paste("PC1 (",varPC1,"%)",sep=""))
p <- p + ylab(paste("PC2 (",varPC2,"%)",sep=""))
p <- p + xlim(-xRange,xRange)
p <- p + ylim(-yRange,yRange)

#p <- p + geom_text_repel(size=3)
#p <- p + theme_classic()
p <- p + theme(axis.text = element_text(size = 7),
               #axis.title = element_text(size = 7, face="bold"),
               axis.title.x = element_text(size=7),
               axis.title.y = element_text(size=7),
               text = element_text(size=7),
               #axis.text.x = element_text(angle=45,hjust=1), 
               #axis.title.x = element_blank(),
               #axis.title.y = element_blank(),
               #panel.border = element_rect(linetype = "solid", colour = "black"),
               panel.background = element_blank(),
               panel.border=element_rect(colour = "black", fill=NA, size=1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),
               axis.ticks=element_line(colour="black"),
               #plot.margin=unit(c(1,1,1,1),"line"),
               plot.margin= margin(10, 10, 5, 5, "pt"),
               #axis.line=element_line(colour="black"),
               plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
               plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
               legend.position="none")

print(p)

#####

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/PCA_plot/03_02_2021"
dir.create(filePath, showWarnings = TRUE, recursive = TRUE)

fileName<-"HCC_PCA_color_by_group_HCC_tumor_vs_HCC_normal_03_02_2021.pdf"

fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=3)
print(p)
dev.off()

######

p <- ggplot(samppca,aes(x=PC1, y=PC2, color=group))
p <- p + scale_color_manual(values=c("steelblue1","tomato","mediumpurple1","grey55"))
#p <- p + geom_point(aes(shape=mutation))
p <- p + geom_hline(yintercept = 0, linetype="dashed",color="grey85")
p <- p + geom_vline(xintercept = 0, linetype="dashed",color="grey85")
p <- p + geom_point(aes(shape=tissue),size=1.5)
p <- p + xlab(paste("PC1 (",varPC1,"%)",sep=""))
p <- p + ylab(paste("PC2 (",varPC2,"%)",sep=""))
p <- p + xlim(-xRange,xRange)
p <- p + ylim(-yRange,yRange)

#p <- p + geom_text_repel(size=3)
#p <- p + theme_classic()
p <- p + theme(axis.text = element_text(size = 7),
               #axis.title = element_text(size = 7, face="bold"),
               axis.title.x = element_text(size=7),
               axis.title.y = element_text(size=7),
               text = element_text(size=7),
               #axis.text.x = element_text(angle=45,hjust=1), 
               #axis.title.x = element_blank(),
               #axis.title.y = element_blank(),
               #panel.border = element_rect(linetype = "solid", colour = "black"),
               panel.background = element_blank(),
               panel.border=element_rect(colour = "black", fill=NA, size=1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),
               axis.ticks=element_line(colour="black"),
               #plot.margin=unit(c(1,1,1,1),"line"),
               plot.margin= margin(10, 10, 5, 5, "pt"),
               #axis.line=element_line(colour="black"),
               plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
               plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
               legend.position="right")

print(p)

#####

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/PCA_plot/03_02_2021"
dir.create(filePath, showWarnings = TRUE, recursive = TRUE)

fileName<-"HCC_PCA_color_by_group_HCC_tumor_vs_HCC_normal_legend_03_02_2021.pdf"

fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=5)
print(p)
dev.off()

