library(ggplot2)
library(ggrepel)
library(plyr)

source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/diffMetaboliteAbundance.R')
source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlot_v3.R')


####
# load lipid name
####
filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/hcc_cancercell"
fileName<-"lipid_types.txt"
#fileName<-"lipid_types_v2.txt"
fileName<-file.path(filePath,fileName)

lipids_definition<-read.table(fileName,sep="\t",header=TRUE,stringsAsFactors = FALSE)

####
# load protected additional clinical information
####

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/PHI"
fileName<-"clin.Rda"
fileName<-file.path(filePath,fileName)
load(fileName)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/PHI"
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

clin_mtDNA_cmoplex_1_mut<-clingen[grepl("complex 1",clingen$Complex..mutation),]$SampleID


#####

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/metabotools_HCC"
#fileName<-"HCC_met_prefiltered_PQN_normalized.Rda"
fileName<-"HCC_met_Elisa_PNQ_normalized_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

# met 737 metabolite x 75 samples
# dim(met)
load(fileName)

#####
# load metabolite info
#####

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_metabolite_info_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

load(file=fileName)

metabolite_lipid<-metinfo[metinfo$SUPER_PATHWAY %in% "Lipid",]


#####
# load smaple info
#####
filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/metabotools_HCC"
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
# HCC, 1            
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
#normalSampleName_hcc<-normalSampleName_hcc[!( normalSampleName_hcc %in% removedSamples)]


#####
# load consensus clustering assignment
#####

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/02_21_2021/combine_mRNA_met_cluster"

fileName<-"met_mRNA_consensus_cluster_assignment.txt"
fileName<-file.path(filePath,fileName)

met_mRNA_cluster_label<-read.table(file=fileName,sep="\t",header=TRUE,stringsAsFactors = FALSE)



#####

samples<-c(tumorSampleName_hcc,normalSampleName_hcc)

met_subset<-met[rownames(met) %in% lipids_definition$name,colnames(met) %in% samples]
#met_subset<-met[rownames(met) %in% metabolite_lipid$name,colnames(met) %in% samples]

met_subset<-met_subset[lipids_definition$name,samples]

######

newMet<-met_subset

rowSD<-apply(newMet,1,sd)

#sdCutoff<-quantile(sort(rowSD,decreasing = TRUE),percentileThreshold)
sdCutoff<-0
rowSubset<-names(rowSD[rowSD>sdCutoff])

numberOfRowsLeft<-length(rowSubset)
numberOfRowsTotal<-length(rownames(newMet))

message(sprintf("Use %s / %s metabolites in the data",numberOfRowsLeft,numberOfRowsTotal))

newMetPCA<-newMet[rowSubset,]
dim(newMetPCA)

met_subset<-newMetPCA

#log2_met_tumor<-log2(met_tumor)

######
filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/02_21_2021/combine_mRNA_met_cluster"

fileName<-"met_mRNA_consensus_cluster_assignment.txt"
fileName<-file.path(filePath,fileName)

met_mRNA_cluster_label<-read.table(file=fileName,sep="\t",header=TRUE,stringsAsFactors = FALSE)
rownames(met_mRNA_cluster_label)<-met_mRNA_cluster_label$sampleName


consensus_sample_order<-rownames(met_mRNA_cluster_label[order(met_mRNA_cluster_label$class,
                                                              met_mRNA_cluster_label$silhouette),])
met_mRNA_cluster_label_ordered<-met_mRNA_cluster_label[consensus_sample_order,]

sample_class_1<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==1,])
sample_class_2<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==2,])
sample_class_3<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==3,])
sample_class_4<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==4,])


######

if(FALSE){

library(survival)
library(survminer)

clingen$sampleName<-rownames(clingen)

cc<-merge(met_mRNA_cluster_label_ordered,clingen,by="sampleName",all.x=TRUE,all.y=TRUE)

mean(cc[cc$class==1,]$PFS_months, na.rm=TRUE)
mean(cc[cc$class==2,]$PFS_months, na.rm=TRUE)
mean(cc[cc$class==3,]$PFS_months, na.rm=TRUE)
mean(cc[cc$class==4,]$PFS_months, na.rm=TRUE)

cc$status<-ifelse(cc$PFS.0.progress..1.no.==0,1,0)

cc<-cc[!is.na(cc$class),]

#fit<-survfit(Surv(PFS_months,PFS.0.progress..1.no.) ~ class, data=cc)
fit<-survfit(Surv(PFS_months,status) ~ WCDChr7, data=cc)
fit<-survfit(Surv(PFS_months,status) ~ LOHorUPD, data=cc)
fit<-survfit(Surv(PFS_months,status) ~ Overall.Stage, data=cc)
fit<-survfit(Surv(PFS_months,status) ~ histology, data=cc)
fit<-survfit(Surv(PFS_months,status) ~ Complex..mutation, data=cc)
fit<-survfit(Surv(PFS_months,status) ~ Somatic.mt.mutation, data=cc)
fit<-survfit(Surv(PFS_months,status) ~ mtDNA_complex, data=cc)


print(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
           )
           #palette = c("#E7B800", "#2E9FDF"))

}

######
# standard error of mean
#sem <- function(x) sd(x)/sqrt(length(x))

output<-list()
test_result<-list()
plotData<-list()

for(idx in 1:nrow(met_subset)){
    
    # 3: palmitate (16:0)
    # 5: stearate (18:0)
    # 9: palmitoleate (16:1n7)
    #11: oleate/vaccenate (18:1)
  
  #https://stackoverflow.com/questions/36691732/how-to-take-log2-of-a-matrix-having-negative-values-for-boxplot-in-r
    abs_log <- function(x){
      x[x==0] <- 1
      si <- sign(x)
      si * log2(si*x)
    }
  
  
    lipid_name<-rownames(met_subset)[idx]
    tumor_vector<-met_subset[idx,tumorSampleName_hcc]
    normal_vector<-met_subset[idx,normalSampleName_hcc]
    
    tumor_vector_class_1<-met_subset[idx,sample_class_1]
    tumor_vector_class_2<-met_subset[idx,sample_class_2]
    tumor_vector_class_3<-met_subset[idx,sample_class_3]
    tumor_vector_class_4<-met_subset[idx,sample_class_4]
    
    
    normal_median<-median(normal_vector)
    #normal_mean<-mean(normal_vector)
    #tumor_vector<-tumor_vector - normal_median
    #normal_vector<-normal_vector - normal_median
    
    log2FC_tumor_vs_normal<-log2(mean(tumor_vector)/mean(normal_vector))
    
    log2FC_tumor_class_1_vs_normal<-log2(mean(tumor_vector_class_1)/mean(normal_vector))
    log2FC_tumor_class_2_vs_normal<-log2(mean(tumor_vector_class_2)/mean(normal_vector))
    log2FC_tumor_class_3_vs_normal<-log2(mean(tumor_vector_class_3)/mean(normal_vector))
    log2FC_tumor_class_4_vs_normal<-log2(mean(tumor_vector_class_4)/mean(normal_vector))
    
    
    
    tumor_vector_log2<-abs_log(tumor_vector)
    normal_vector_log2<-abs_log(normal_vector)
    
    tumor_vector_relative_to_normal<-(tumor_vector)/(mean(normal_vector))
    normal_vector_relative_to_normal<-(normal_vector)/(mean(normal_vector))
    #normal_vector<-abs_log(normal_vector)
    
    tumor_mean<-mean(tumor_vector)
    tumor_sd<-sd(tumor_vector)
    tumor_sem<-sd(tumor_vector)/sqrt(length(tumor_vector))
    
    normal_mean<-mean(normal_vector)
    normal_sd<-sd(normal_vector)
    normal_sem<-sd(normal_vector)/sqrt(length(normal_vector))
    
    res.wilcox<-wilcox.test(tumor_vector,normal_vector)
    wilcox_pvalue<-res.wilcox$p.value
    
    #log2FC<-mean(tumor_vector)/mean(normal_vector)
    
    res.t<-t.test(tumor_vector,normal_vector)
    t_test_pvalue<-res.t$p.value
    
    group<-c("tumor","normal")
    lipid<-rep(lipid_name,2)
    mean<-c(tumor_mean,normal_mean)
    sd<-c(tumor_sd,normal_sd)
    sem<-c(tumor_sem,normal_sem)
    
    output[[idx]]<-data.frame(group,lipid,mean,sem,stringsAsFactors = FALSE)
    
    if(FALSE){
    plotData[[idx]]<-data.frame(lipid=rep(lipid_name,length(c(tumor_vector_log2,normal_vector_log2))),
                                value=c(tumor_vector_log2,normal_vector_log2),
                                group=c(rep("tumor",length(tumor_vector_log2)),rep("normal",length(normal_vector_log2))),
                                stringsAsFactors = FALSE)
    }
    
    plotData[[idx]]<-data.frame(lipid=rep(lipid_name,length(c(tumor_vector_relative_to_normal,normal_vector_relative_to_normal))),
                                value=c(tumor_vector_relative_to_normal,normal_vector_relative_to_normal),
                                group=c(rep("tumor",length(tumor_vector_relative_to_normal)),rep("normal",length(normal_vector_relative_to_normal))),
                                stringsAsFactors = FALSE)
    
    test_result[[idx]]<-data.frame(lipid_name,
                                   log2FC_tumor_vs_normal,
                                   log2FC_tumor_class_1_vs_normal,
                                   log2FC_tumor_class_2_vs_normal,
                                   log2FC_tumor_class_3_vs_normal,
                                   log2FC_tumor_class_4_vs_normal,
                                   wilcox_pvalue,t_test_pvalue,stringsAsFactors = FALSE)
    
}

output<-rbind.fill(output)

plotData<-rbind.fill(plotData)

test_result<-rbind.fill(test_result)
rownames(test_result)<-test_result$lipid_name

test_result$padj_wilcox<-p.adjust(test_result$wilcox_pvalue,method="BH")
test_result$log2FC<-test_result$log2FC_tumor_vs_normal

plot<-makeVolcanoPlot(res=test_result,pvalueCol = 'padj_wilcox',threshold_pvalue = 0.1,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/03_16_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_vs_normal_log2_scale_volcanoplot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()



#library(ggpubr)
#res<-compare_means(mean~group,data=output)


abs_log <- function(x){
  x[x==0] <- 1
  si <- sign(x)
  si * log2(si*x)
}

output$log2_mean<-abs_log(output$mean)


#####
# volcano plot
#####
# work further here


test_result$group2<-NA
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "SFA",]$name,]$group2<-"SFA"
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "MUFA",]$name,]$group2<-"MUFA"
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "PUFA",]$name,]$group2<-"PUFA"


threshold_log2fc<-1
threshold_pvalue<-0.05
pvalueCol<-"padj_wilcox"

res<-test_result
res$log10pvalue<-(-log10(res$padj_wilcox))
res$label<-res$lipid_name
foldChangeRange<-ceiling(max(abs(res$log2FC)))

foldChangeBreaks<-seq(-foldChangeRange,foldChangeRange,1)

yAxisRange<-ceiling(max(abs(res$log10pvalue)))+1
yAxisBreaks<-seq(0,yAxisRange,1)



p <- ggplot(res,aes(x=log2FC,y=log10pvalue,label = label)) + geom_point(aes(color=group2))
#p <- p + theme_minimal(base_size = 16) 
#p <- p + xlim(-5,5)
p <- p + geom_text_repel(
  #data=res[abs(res$log2FC)>threshold_log2fc & res[,pvalueCol] < threshold_pvalue,],
  data=res,
  size = 2,
  min.segment.length = 0,
  box.padding = 0,
  max.overlaps = Inf,
  nudge_x = 0.15,
  #nudge_y = 0.15,
  #segment.curvature = -0.1,
  #segment.ncp = 3,
  #segment.angle = 20
)   
p <- p + xlab('Log2 Fold Change') 
p <- p + ylab('-log10 (p-value)') 
p <- p + geom_rug(alpha=1/2, position="jitter")
#p <- p + scale_color_manual(name = 'Differential Abundance',values = c('Decrease' = 'blue','Increase' = 'red','None' = 'gray','Black' = 'Black'))
#p <- p + theme(legend.position = 'bottom')
p <- p + geom_hline(yintercept = -log10(threshold_pvalue), linetype="dashed", color="gray")
p <- p + geom_vline(xintercept = abs(threshold_log2fc), linetype="dashed", color="gray")
p <- p + geom_vline(xintercept = -abs(threshold_log2fc), linetype="dashed", color="gray")
#p <- p + xlim(-foldChangeRange,foldChangeRange)
p <- p + scale_x_continuous(breaks=foldChangeBreaks)
#p <- p + ylim(0,max(res$log10pvalue)+5)
p <- p + expand_limits(x=c(-foldChangeRange,foldChangeRange),y=c(0,yAxisRange))
#p <- p + scale_y_continuous(expand = expansion(mult=c(0, 0.1),add=1))
p <- p + scale_y_continuous(#expand = c(0,0),
  breaks=yAxisBreaks)
p <- p + theme_classic()
p <- p + theme(axis.text = element_text(size = 7, face="bold"),
               #axis.title = element_text(size = 7, face="bold"),
               axis.title.x = element_text(size=7),
               axis.title.y = element_text(size=7),
               text = element_text(size=7),
               #axis.text.x = element_text(angle=45,hjust=1), 
               #axis.title.x = element_blank(),
               #axis.title.y = element_blank(),
               #panel.border = element_rect(linetype = "solid", colour = "black"),
               panel.border = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line=element_line(colour="black"),
               plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
               plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
               legend.position="bottom")
print(p)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/03_16_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_vs_normal_log2_scale_color_volcanoplot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=4,width=4)
print(p)
dev.off()

#####
# make log2FC plot
#####

test_result$group2<-NA
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "SFA",]$name,]$group2<-"SFA"
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "MUFA",]$name,]$group2<-"MUFA"
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "PUFA",]$name,]$group2<-"PUFA"


res<-test_result
res$log10pvalue<-(-log10(res$padj_wilcox))
res$label<-res$lipid_name
foldChangeRange<-ceiling(max(abs(res$log2FC)))

foldChangeBreaks<-seq(-foldChangeRange,foldChangeRange,1)

yAxisRange<-ceiling(max(abs(res$log10pvalue)))+1
yAxisBreaks<-seq(0,yAxisRange,1)

res$lipid_name<-factor(res$lipid_name,levels=lipids_definition$name)

graph <- ggplot(res,aes(x=lipid_name,y=log2FC, fill=group2))
graph <- graph + geom_bar(stat="identity")
graph <- graph + xlab("Lipids")
graph <- graph + ylab("log2 fold change (tumor/normal)")
graph <- graph + scale_fill_manual(name = 'Lipid',values = c('SFA' = 'grey','MUFA' = 'royalblue','PUFA' = 'red'))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.major.y = element_blank(), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "bottom"
)

print(graph)

######
# make relative fold change
#####

significance_data<-test_result[,c("lipid_name","log2FC_tumor_vs_normal","wilcox_pvalue")]
significance_data$x_coord<-seq(1,nrow(significance_data),1)
significance_data$y_coord<-2^(significance_data$log2FC_tumor_vs_normal)+0.8

significance_data$label<-""
for(idx in 1:nrow(significance_data)){
  if(significance_data$wilcox_pvalue[idx]<0.05){significance_data$label[idx]<-"*"}
  if(significance_data$wilcox_pvalue[idx]<0.01){significance_data$label[idx]<-"**"}
  if(significance_data$wilcox_pvalue[idx]<0.001){significance_data$label[idx]<-"***"}
}


foldChangeRange<-ceiling(max(abs(plotData$value)))
foldChangeRange<-15
foldChangeBreaks<-seq(-foldChangeRange,foldChangeRange,1)

#yAxisRange<-ceiling(max(abs(res$log10pvalue)))+1
yAxisRange<-foldChangeRange
yAxisBreaks<-seq(0,yAxisRange,1)

plotData$group2<-NA
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "SFA",]$name,]$group2<-"SFA"
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "MUFA",]$name,]$group2<-"MUFA"
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "PUFA",]$name,]$group2<-"PUFA"


plotData$lipid<-factor(plotData$lipid,levels=lipids_definition$name)
#plotData$lipid_label<-lipids_definition$label
#plotData$lipid_label<-factor(plotData$lipid_label,levels=lipids_definition$label)

colnames(plotData)[1]<-"name"
plotData2<-merge(plotData,lipids_definition,by="name",all.x=TRUE)
plotData2$label<-factor(plotData2$label,levels=lipids_definition$label)


graph <- ggplot(plotData2,aes(x=label,y=value,fill=group))
#graph <- graph + geom_bar(stat="identity")
graph <- graph + stat_summary(aes(y=value),
                              fun="mean",geom="bar",
                              position=position_dodge(width=0.7),
                              stat="identity", 
                              width=0.6)
graph <- graph + stat_summary(fun.data = "mean_se", 
                              geom = "errorbar",
                              position=position_dodge(width=0.7),
                              stat="identity",
                              size=0.3,
                              width=0.2
                              )

graph <- graph + geom_text(data=significance_data, aes(x=x_coord, y=y_coord,label=label,fill=NULL,angle=90,vjust=0.7), size=(7 / .pt))
graph <- graph + xlab("Lipids")
graph <- graph + ylab("Relative ratio")
#graph <- graph + scale_fill_manual(name = 'Lipid',values = c('SFA' = 'grey','MUFA' = 'royalblue','PUFA' = 'red'))
graph <- graph + scale_fill_manual(name = 'Lipid',values = c('normal' = 'grey','tumor' = 'indianred1'))
graph <- graph + scale_y_continuous(expand = expansion(mult=c(0, 0.1)),
                             breaks=yAxisBreaks)
graph <- graph + theme_classic()
graph <- graph + theme(axis.text = element_text(size = 7),
               #axis.title = element_text(size = 7, face="bold"),
               axis.title.x = element_text(size=7),
               axis.title.y = element_text(size=7),
               text = element_text(size=7),
               axis.text.x=element_text(colour="black", angle=90, hjust=1,vjust=0.5),
               axis.text.y=element_text(colour="black"),
               axis.ticks=element_line(colour="black"),
               #axis.text.x = element_text(angle=45,hjust=1), 
               #axis.title.x = element_blank(),
               #axis.title.y = element_blank(),
               #panel.border = element_rect(linetype = "solid", colour = "black"),
               panel.border = element_blank(),
               panel.grid.major = element_line(colour = "grey",size=0.1), 
               #panel.grid.major = element_blank(), 
               #panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line=element_line(colour="black"),
               plot.margin= margin(t = 10, r = 5, b = 5, l = 10, "pt"),
               plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
               plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
               legend.text=element_text(size=7),
               legend.key.size = unit(3,"mm"),
               legend.position="bottom")

print(graph)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/03_16_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_vs_normal_relativeRatio_barplot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=4)
print(graph)
dev.off()

######




#####
# make log2FC plot
#####
test_result$group2<-NA
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "SFA",]$name,]$group2<-"SFA"
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "MUFA",]$name,]$group2<-"MUFA"
test_result[test_result$lipid_name %in% lipids_definition[lipids_definition$type == "PUFA",]$name,]$group2<-"PUFA"


res<-test_result

d0<-data.frame(test_result[,c("lipid_name","log2FC_tumor_vs_normal")],rep("tumor",nrow(test_result)),stringsAsFactors = FALSE)
colnames(d0)<-c("lipid","log2FC","type")
d1<-data.frame(test_result[,c("lipid_name","log2FC_tumor_class_1_vs_normal")],rep("tumor_class_1",nrow(test_result)),stringsAsFactors = FALSE)
colnames(d1)<-c("lipid","log2FC","type")
d2<-data.frame(test_result[,c("lipid_name","log2FC_tumor_class_2_vs_normal")],rep("tumor_class_2",nrow(test_result)),stringsAsFactors = FALSE)
colnames(d2)<-c("lipid","log2FC","type")
d3<-data.frame(test_result[,c("lipid_name","log2FC_tumor_class_3_vs_normal")],rep("tumor_class_3",nrow(test_result)),stringsAsFactors = FALSE)
colnames(d3)<-c("lipid","log2FC","type")
d4<-data.frame(test_result[,c("lipid_name","log2FC_tumor_class_4_vs_normal")],rep("tumor_class_4",nrow(test_result)),stringsAsFactors = FALSE)
colnames(d4)<-c("lipid","log2FC","type")



  
res<-rbind( d0,d1,d2,d3,d4
)

res$group2<-NA
res[res$lipid %in% lipids_definition[lipids_definition$type == "SFA",]$name,]$group2<-"SFA"
res[res$lipid %in% lipids_definition[lipids_definition$type == "MUFA",]$name,]$group2<-"MUFA"
res[res$lipid %in% lipids_definition[lipids_definition$type == "PUFA",]$name,]$group2<-"PUFA"

res$lipid<-factor(res$lipid,levels=lipids_definition$name)

graph <- ggplot(res,aes(x=lipid,y=log2FC, fill=type))
graph <- graph + geom_bar(stat="identity",position=position_dodge(),color="black",width=0.5)
graph <- graph + xlab("Lipids")
graph <- graph + ylab("log2 fold change (tumor/normal)")
#graph <- graph + scale_fill_manual(name = 'Lipid',values = c('SFA' = 'grey','MUFA' = 'royalblue','PUFA' = 'red'))
graph <- graph + scale_fill_manual(name = 'Lipid',values = c('tumor' = 'indianred2',
                                                             'tumor_class_1' = 'red',
                                                             'tumor_class_2' = 'royalblue2',
                                                             'tumor_class_3' = 'mediumpurple1',
                                                             'tumor_class_4' = 'pink1'))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.major.y = element_blank(), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "bottom"
)

print(graph)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/02_23_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_classes_log2fc_barplot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=7,width=15)
print(graph)
dev.off()

######

graph <- ggplot(res,aes(x=lipid_name,y=log2FC, fill=group2))
graph <- graph + geom_bar(stat="identity")
graph <- graph + geom_rug()
graph <- graph + xlab("Lipids")
graph <- graph + ylab("log2 fold change (tumor/normal)")
graph <- graph + scale_fill_manual(name = 'Lipid',values = c('SFA' = 'grey','MUFA' = 'royalblue','PUFA' = 'red'))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.major.y = element_blank(), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "bottom"
)

print(graph)




######


filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/02_23_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_vs_normal_log2fc_barplot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=7,width=15)
print(graph)
dev.off()



#####
# original boxplot in log2(metabolite abundance) scale
#####

if(TRUE){

plotData$lipid<-factor(plotData$lipid,levels=lipids_definition$name)

plotData$group2<-NA
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "SFA",]$name,]$group2<-"SFA"
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "MUFA",]$name,]$group2<-"MUFA"
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "PUFA",]$name,]$group2<-"PUFA"

yRangeMin<-round(min(plotData$value)-1)
yRangeMax<-round(max(plotData$value)+1)

graph<-ggplot(data=plotData,aes(x=lipid,y=value,fill=group))
graph <- graph + geom_boxplot(outlier.shape = NA)
#graph <- graph + geom_beeswarm()
graph <- graph + xlab("Lipids")
graph <- graph + ylab("log2(Metabolite abundance)")
graph <- graph + scale_fill_manual(values=c("royalblue","red"))
graph <- graph + expand_limits(y=c(yRangeMin,yRangeMax))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "grey",size=0.1), 
                       panel.grid.major.y = element_blank(), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/03_16_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_vs_normal_log2_scale_boxplot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=7,width=15)
print(graph)
dev.off()



}




#####

plotData$lipid<-factor(plotData$lipid,levels=lipids_definition$name)

plotData$group2<-NA
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "SFA",]$name,]$group2<-"SFA"
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "MUFA",]$name,]$group2<-"MUFA"
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "PUFA",]$name,]$group2<-"PUFA"

yRangeMin<-round(min(plotData$value)-1)
yRangeMax<-round(max(plotData$value)+1)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

graph<-ggplot(data=plotData,aes(x=lipid,y=value,fill=group))
#graph <- graph + geom_boxplot(outlier.shape = NA)
graph <- graph + geom_violin()
graph <- graph + stat_summary(fun.data=data_summary)
#graph <- graph + facet_wrap(~group2, ncol=1)
#graph <- graph + geom_beeswarm()
graph <- graph + xlab("Lipids")
graph <- graph + ylab("log2(Metabolite abundance)")
graph <- graph + scale_fill_manual(values=c("royalblue","red"))
graph <- graph + expand_limits(y=c(yRangeMin,yRangeMax))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "grey",size=0.1), 
                       panel.grid.major.y = element_blank(), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/02_23_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_vs_normal_log2_scale_boxplot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=7,width=15)
print(graph)
dev.off()

######

library(cowplot)


plotData$lipid<-factor(plotData$lipid,levels=lipids_definition$name)

plotData$group2<-NA
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "SFA",]$name,]$group2<-"SFA"
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "MUFA",]$name,]$group2<-"MUFA"
plotData[plotData$lipid %in% lipids_definition[lipids_definition$type == "PUFA",]$name,]$group2<-"PUFA"

yRangeMin<-round(min(plotData$value)-1)
yRangeMax<-round(max(plotData$value)+1)

plotData_SFA<-plotData[plotData$group2 %in% "SFA",]

graph<-ggplot(data=plotData_SFA,aes(x=lipid,y=value,fill=group))
graph <- graph + geom_boxplot(outlier.shape = NA)
graph <- graph + facet_wrap(~group2, ncol=1)
#graph <- graph + geom_beeswarm()
graph <- graph + xlab("Lipids")
graph <- graph + ylab("log2(Metabolite abundance)")
graph <- graph + scale_fill_manual(values=c("royalblue","red"))
graph <- graph + expand_limits(y=c(yRangeMin,yRangeMax))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)
p1<-graph

####

plotData_MUFA<-plotData[plotData$group2 %in% "MUFA",]

graph<-ggplot(data=plotData_MUFA,aes(x=lipid,y=value,fill=group))
graph <- graph + geom_boxplot(outlier.shape = NA)
graph <- graph + facet_wrap(~group2, ncol=1)
#graph <- graph + geom_beeswarm()
graph <- graph + xlab("Lipids")
graph <- graph + ylab("log2(Metabolite abundance)")
graph <- graph + scale_fill_manual(values=c("royalblue","red"))
graph <- graph + expand_limits(y=c(yRangeMin,yRangeMax))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)
p2<-graph


####

plotData_PUFA<-plotData[plotData$group2 %in% "PUFA",]

graph<-ggplot(data=plotData_PUFA,aes(x=lipid,y=value,fill=group))
graph <- graph + geom_boxplot(outlier.shape = NA)
graph <- graph + facet_wrap(~group2, ncol=1)
#graph <- graph + geom_beeswarm()
graph <- graph + xlab("Lipids")
graph <- graph + ylab("log2(Metabolite abundance)")
graph <- graph + scale_fill_manual(values=c("royalblue","red"))
graph <- graph + expand_limits(y=c(yRangeMin,yRangeMax))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)
p3<-graph

#####

plot_row<-plot_grid(p1,p2,p3)

plot_grid(plot_row,ncol=1)

plot_grid(p1,p2,p3, nrow=1)


print(graph)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/02_23_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_vs_normal_log2_scale.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=7,width=15)
print(graph)
dev.off()


#####
# add significance "*" mark to barplot
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
#
output$lipid<-factor(output$lipid,levels=lipids_definition$name)

graph <- ggplot(data=output,aes(x=lipid,y=log2_mean,fill=group))
graph <- graph + geom_bar(stat="identity", position=position_dodge())
#graph <- graph + geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), 
#                                 width=.2,
#                                 position=position_dodge(.9))

graph <- graph + xlab("Lipids")
graph <- graph + ylab("Metabolite abundance")
#graph <- graph + scale_fill_brewer(palette="Paired")
graph <- graph + scale_fill_manual(values=c("royalblue","red"))
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=12, color="black"),
                       axis.text = element_text(size = 12, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.text.x = element_text(size=12, angle=45, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       #panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)


filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/lipid_metabolism/02_22_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_lipid_abundance_SFA_MUFA_PUFA_without_acyl_carnitine_tumor_vs_normal_log2_scale.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=7,width=15)
print(graph)
dev.off()





