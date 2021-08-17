library(ggplot2)
library(ggrepel)

source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/diffMetaboliteAbundance.R')
source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlot_v3.R')
source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlot_for_publication.R')


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

clin_mtDNA_cmoplex_1_mut<-clingen[grepl("complex 1",clingen$Complex..mutation),]$SampleID


#####

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


# 7 samples in the PD and PTCTV with adjacent normal samples
pd_and_ptctv_normal_group<-c("NT_PDTC_Mut","NT_PDTC_WT","NT_PTC_WT")
normalSampleName_pd_and_ptctv<-sampleInfo_Normal[sampleInfo_Normal$Group %in% pd_and_ptctv_normal_group,]$CLIENT_IDENTIFIER


# 4 samples in the HA with adjacent normal samples
ha_normal_group<-c("NT_HA")
normalSampleName_ha<-sampleInfo_Normal[sampleInfo_Normal$Group %in% ha_normal_group,]$CLIENT_IDENTIFIER


######

matchedNormalSampleName_hcc<-sampleInfo_Normal[sampleInfo_Normal$Group %in% hcc_normal_group,]$CLIENT_IDENTIFIER
matchedTumorSampleName_hcc<-sampleInfo_Normal[sampleInfo_Normal$Group %in% hcc_normal_group,]$MATCHED

# drop YD29 (HCC HMIN), YD65 (PDTC) tumor samples at the beginning
# drop YD30, YD66 matched normal samples for analysis    
removedSamples<-c("YD29","YD65","YD30","YD66")
matchedNormalSampleName_hcc<-matchedNormalSampleName_hcc[!( matchedNormalSampleName_hcc %in% removedSamples)]
matchedTumorSampleName_hcc<-matchedTumorSampleName_hcc[!( matchedTumorSampleName_hcc %in% removedSamples)]
#normalSampleName_hcc<-normalSampleName_hcc[!( normalSampleName_hcc %in% removedSamples)]




#####

if(FALSE){
library(limma)

edata<-log2(met)


sampleInfo<-sampleInfo[sampleInfo$CLIENT_IDENTIFIER %in% colnames(edata),]
sampleInfo$TISSUE<-factor(sampleInfo$TISSUE,levels=c("NORMAL","TUMOR"))

model<-model.matrix(~0+TISSUE,data=sampleInfo)
colnames(model)<-c("normal","tumor")

fit<-lmFit(edata,model)

contrast.matrix<-makeContrasts(tumorVSnormal = tumor - normal, levels=model)

fit2<-contrasts.fit(fit, contrast.matrix)

fit2<-eBayes(fit2)

results_tumor_vs_normal_limma_test<-topTable(fit2, number=nrow(fit2), adjust.method = "BH")

results_tumor_vs_normal_limma_test$metabolite<-rownames(results_tumor_vs_normal_limma_test)
results_tumor_vs_normal_wilcox_test$metabolite<-rownames(results_tumor_vs_normal_wilcox_test)

cc<-merge(results_tumor_vs_normal_limma_test,results_tumor_vs_normal_wilcox_test,by="metabolite")

}

#####
#source('~/work/Ed_lab/code_dev/metabolomics_pipeline/HCC/diffMetaboliteAbundance.R')
#source('~/work/Ed_lab/code_dev/metabolomics_pipeline/HCC/makeVolcanoPlotPublication.R')




if(FALSE){
# log2( mean(ix2) / mean (ix1) ) fold change
results_tumor_vs_normal_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=normalSampleName,ix2=tumorSampleName)
results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

head(results_tumor_vs_normal_wilcox_test)

####
res<-results_tumor_vs_normal_wilcox_test

plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.1,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/result/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_prefiltered_PQN_normalized_tumor_vs_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/result/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_vs_normal.csv"
fileName<-file.path(filePath,fileName)

res$name<-rownames(res)

write.table(res,file=fileName,sep=",",quote=FALSE,row.names=TRUE,col.names = TRUE)



####
res<-result_tumor_vs_normal
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']


# minor correction
metinfo[metinfo$H_SUPER_PATHWAY %in% "Amino acid",]$H_SUPER_PATHWAY<-"Amino Acid"
metinfo[metinfo$H_SUB_PATHWAY %in% "Fatty acid metabolism (also BCAA metabolism)",]$H_SUB_PATHWAY<-"Fatty Acid Metabolism (also BCAA Metabolism)"
metinfo[metinfo$H_SUB_PATHWAY %in% "Vitamin B6 metabolism",]$H_SUB_PATHWAY<-"Vitamin B6 Metabolism"
metinfo[metinfo$H_SUB_PATHWAY %in% "Riboflavin metabolism",]$H_SUB_PATHWAY<-"Riboflavin Metabolism"
metinfo[metinfo$H_SUB_PATHWAY %in% "Polyamine metabolism",]$H_SUB_PATHWAY<-"Polyamine Metabolism"



res$H_COMP_IDstr<-metinfo[rownames(res),'H_COMP_IDstr']
res$H_SUPER_PATHWAY<-metinfo[rownames(res),'H_SUPER_PATHWAY']
res$H_SUB_PATHWAY<-metinfo[rownames(res),'H_SUB_PATHWAY']

###


#super_pathway_list<-unique(metinfo$H_SUPER_PATHWAY)

# only pick 7 super pathway annotated from metabolon
super_pathway_list<-c("Lipid",
                      "Cofactors and Vitamins",
                      "Amino Acid",
                      "Carbohydrate",
                      "Energy",
                      "Nucleotide",
                      "Peptide")

#test_pathway<-c("Peptide","Energy")




tmpFrame<-res
#tmpFrame<-res[res$H_SUPER_PATHWAY %in% test_pathway,]
tmpFrame<-res[res$H_SUPER_PATHWAY %in% super_pathway_list,]
tmpFrame<-tmpFrame[!is.na(tmpFrame$H_SUB_PATHWAY),]

#sub_pathway_list<-unique(tmpFrame$H_SUB_PATHWAY)
order_list<-list()

for(idx in 1:length(super_pathway_list)){
  tmpVector<-unique(tmpFrame[tmpFrame$H_SUPER_PATHWAY %in% super_pathway_list[idx],]$H_SUB_PATHWAY)
  order_list[[idx]]<-tmpVector
}

order_list<-unlist(order_list)

super_pathway_frame<-list()

for(idx in 1:length(super_pathway_list)){
  
     #idx<-7
     #message(sprintf("idx1: %s",idx)) 
     currentPathway<-super_pathway_list[idx]  
     step1<-tmpFrame[tmpFrame$pathway %in% currentPathway,]
     
     sub_pathway_list<-unique(step1$H_SUB_PATHWAY)
     sub_pathway_frame<-list()
     
     for( idx2 in 1:length(sub_pathway_list)){
        
         #message(sprintf("idx2: %s",idx2))  
         #idx2<-1
         currentSubPathway<-sub_pathway_list[idx2]
         step2<-step1[step1$H_SUB_PATHWAY %in% currentSubPathway,]
         numOfmetabolites<-nrow(step2)
         step2$subPathwayName<-paste(step2$H_SUB_PATHWAY," (",numOfmetabolites,")",sep="")
         
         step2_other<-tmpFrame[!(tmpFrame$pathway %in% currentPathway),]
         log2FC_reference<-mean(step2_other$log2FC)
         step2$log2FC_relative<-step2$log2FC - log2FC_reference
         sub_pathway_frame[[idx2]]<-step2
     }
     
     step1_merged_frame<-rbind.fill(sub_pathway_frame)
     super_pathway_frame[[idx]]<-step1_merged_frame
}

inputFrame<-rbind.fill(super_pathway_frame)

#inputFrame$H_SUB_PATHWAY<-factor(inputFrame$H_SUB_PATHWAY,levels=rev(order_list))

order_list<-list()

for(idx in 1:length(super_pathway_list)){
  tmpVector<-unique(inputFrame[inputFrame$H_SUPER_PATHWAY %in% super_pathway_list[idx],]$subPathwayName)
  order_list[[idx]]<-tmpVector
}

order_list<-unlist(order_list)

inputFrame$subPathwayName<-factor(inputFrame$subPathwayName,levels=rev(order_list))


#inputFrame<-tmpFrame

#graph <- ggplot(inputFrame, aes(x=H_SUB_PATHWAY,y=log2FC_relative,fill=H_SUPER_PATHWAY))
graph <- ggplot(inputFrame, aes(x=subPathwayName,y=log2FC_relative,fill=H_SUPER_PATHWAY))
graph <- graph + geom_boxplot()
#graph <- graph + geom_violin()
graph <- graph + geom_hline(yintercept=0,color="black",size=0.6)
graph <- graph + coord_flip()
graph <- graph + theme_classic()
graph <- graph + theme(axis.text = element_text(size = 7),
                  #axis.title = element_text(size = 7, face="bold"),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size=10),
                  #axis.title.x = element_blank(),
                  #axis.title.y = element_blank(),
                  #panel.border = element_rect(linetype = "solid", colour = "black"),
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line=element_line(colour="black"),
                  plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                  plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                  #legend.position="none")
                  legend.position = "right"
                  )
#legend.position = c(0.65,0.2),
#legend.background = element_rect(color = NA, 
#                                 fill = "transparent", size = 1, linetype = "solid"),
#legend.text = element_text(size = 7, colour = "black", angle = 0))

print(graph)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath, showWarnings = TRUE, recursive = TRUE)
#fileName<-"HCC_prefiltered_PQN_normalized_PCA_by_histology.pdf"
fileName<-"HCC_prefiltered_PQN_normalized_tumor_vs_normal_log2FC_relative_plot.pdf"
#fileName<-"AML_merged_normalized_heatmap_by_group_time_late.pdf"
#fileName<-"AML_merged_normalized_heatmap_by_group_time_late.pdf"

fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(graph)
dev.off()


}

#######
# set0

# log2( mean(ix2) / mean (ix1) ) fold change
results_tumor_vs_normal_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=normalSampleName,ix2=tumorSampleName)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

####
res<-results_tumor_vs_normal_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###
plot<-makeVolcanoPlotPublication(res=res,
                                 pvalueCol = 'padj',
                                 threshold_pvalue = 0.05,
                                 threshold_log2fc = 1,
                                 numOfCandidateLabels = 5,
                                 show_legend_position = "none")

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_tumor_vs_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_vs_normal.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#######
# set1

# log2( mean(ix2) / mean (ix1) ) fold change
result_hcc_vs_pd<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_pd,ix2=tumorSampleName_hcc)

####
res<-result_hcc_vs_pd
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_Elisa_PQN_normalized_hcc_vs_pd_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_vs_pd.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#######
# set2

# log2( mean(ix2) / mean (ix1) ) fold change
result_hcc_vs_ptctv<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_ptctv,ix2=tumorSampleName_hcc)

####
res<-result_hcc_vs_ptctv
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_Elisa_PQN_normalized_hcc_vs_ptctv_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_vs_ptctv.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#######
#set 3
# log2( mean(ix2) / mean (ix1) ) fold change
result_pd_vs_ptctv<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_ptctv,ix2=tumorSampleName_pd)

####
res<-result_pd_vs_ptctv
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_Elisa_PQN_normalized_pd_vs_ptctv_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_pd_vs_ptctv.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#######
# set 4
#######
# log2( mean(ix2) / mean (ix1) ) fold change
result_HCC_vs_pd_and_ptctv<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_pd_and_ptctv,ix2=tumorSampleName_hcc)

####
res<-result_HCC_vs_pd_and_ptctv
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

plot<-makeVolcanoPlotPublication(res=res,
                                 pvalueCol = 'padj',
                                 threshold_pvalue = 0.05,
                                 threshold_log2fc = 1,
                                 numOfCandidateLabels = 5,
                                 show_legend_position = "none")

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_hcc_vs_pd_and_ptctv_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_vs_pd_and_ptctv.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#####
# set 5
#####
# log2( mean(ix2) / mean (ix1) ) fold change
result_hcc_HWIDE_vs_hcc_HMIN<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_hcc_HMIN,ix2=tumorSampleName_hcc_HWIDE)

####
res<-result_hcc_HWIDE_vs_hcc_HMIN
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_hcc_HWIDE_vs_hcc_HMIN_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_HWIDE_vs_hcc_HMIN.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#######
# set6
######

# log2( mean(ix2) / mean (ix1) ) fold change
results_hcc_tumor_vs_hcc_normal_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=normalSampleName_hcc,ix2=tumorSampleName_hcc)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

####
res<-results_hcc_tumor_vs_hcc_normal_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###
plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1, numOfCandidateLabels=20)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_hcc_tumor_vs_hcc_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_tumor_vs_hcc_normal.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)


#######
# set7
######

# log2( mean(ix2) / mean (ix1) ) fold change
results_hcc_vs_ha<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_ha,ix2=tumorSampleName_hcc)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

####
res<-results_hcc_vs_ha
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

selected_labels<-c("pyruvate",
                   "2-phosphoglycerate",
                   "aspartate",
                   "phosphoenolpyruvate (PEP)"
                   #"citrate",
                   #"aconitate [cis or trans]"
                   )

plot<-makeVolcanoPlotPublication(res=res,
                                 pvalueCol = 'padj',
                                 threshold_pvalue = 0.05,
                                 threshold_log2fc = 1,
                                 show_labels = selected_labels,
                                 show_legend_position = "none")

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_Elisa_PQN_normalized_hcc_vs_ha_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=3.5)
print(plot)
dev.off()


plot<-makeVolcanoPlotPublication(res=res,
                                 pvalueCol = 'padj',
                                 threshold_pvalue = 0.05,
                                 threshold_log2fc = 1,
                                 show_labels = selected_labels,
                                 show_legend_position = "bottom")

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_Elisa_PQN_normalized_hcc_vs_ha_volcano_plot_legend.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=3.5)
print(plot)
dev.off()




filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_vs_ha.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)


#######
# set8

# log2( mean(ix2) / mean (ix1) ) fold change
results_matched_hcc_tumor_vs_matched_hcc_normal_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=matchedNormalSampleName_hcc,ix2=matchedTumorSampleName_hcc)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

####
res<-results_matched_hcc_tumor_vs_matched_hcc_normal_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###
plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

#plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_matched_hcc_tumor_vs_matched_hcc_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_matched_hcc_tumor_vs_matched_hcc_normal.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

########

#######
# set9

tumorSampleName_hcc_no_recurrent<-c(tumorSampleName_hcc_HMIN,tumorSampleName_hcc_HWIDE_no_recurrent)

# log2( mean(ix2) / mean (ix1) ) fold change
results_matched_hcc_HWIDE_recurrent_vs_hcc_no_recurrent_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_hcc_no_recurrent,ix2=tumorSampleName_hcc_HWIDE_recurrent)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

#x<-met
#ix1<-tumorSampleName_hcc_no_recurrent
#ix2<-tumorSampleName_hcc_HWIDE_recurrent

####
res<-results_matched_hcc_HWIDE_recurrent_vs_hcc_no_recurrent_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

#res<-res[!(res$pathway %in% "Xenobiotics"),]
#res$padj2<-p.adjust(res$pvalue,method="BH")

#source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlotPublication_v3.R')


plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

#plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_hcc_HWIDE_recurrent_vs_hcc_no_recurrent_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_HWIDE_recurrent_vs_hcc_no_recurrent.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)


######


# log2( mean(ix2) / mean (ix1) ) fold change
results_matched_hcc_HWIDE_recurrent_vs_hcc_normal_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=normalSampleName_hcc,ix2=tumorSampleName_hcc_HWIDE_recurrent)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

#x<-met
#ix1<-tumorSampleName_hcc_no_recurrent
#ix2<-tumorSampleName_hcc_HWIDE_recurrent

####
res<-results_matched_hcc_HWIDE_recurrent_vs_hcc_normal_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

#res<-res[!(res$pathway %in% "Xenobiotics"),]
#res$padj2<-p.adjust(res$pvalue,method="BH")

#source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlotPublication_v3.R')


plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

#plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_hcc_HWIDE_recurrent_vs_hcc_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_HWIDE_recurrent_vs_hcc_normal.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)


######

#######
# set10

# 17 HCC tumors without mtDNA complex 1 mutation
tumorSampleName_hcc_no_mtDNA_complex_1_mut<-setdiff(tumorSampleName_hcc,clin_mtDNA_cmoplex_1_mut)

# 15 HCC tumors with mtDNA complex 1 mutation

# log2( mean(ix2) / mean (ix1) ) fold change
results_matched_hcc_mtDNA_c1_mut_vs_hcc_no_mtDNA_c1_mut_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_hcc_no_mtDNA_complex_1_mut,ix2=clin_mtDNA_cmoplex_1_mut)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

#x<-met
#ix1<-tumorSampleName_hcc_no_recurrent
#ix2<-tumorSampleName_hcc_HWIDE_recurrent

####
res<-results_matched_hcc_mtDNA_c1_mut_vs_hcc_no_mtDNA_c1_mut_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

#res<-res[!(res$pathway %in% "Xenobiotics"),]
#res$padj2<-p.adjust(res$pvalue,method="BH")

#source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlotPublication_v3.R')

plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

#plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_matched_hcc_tumor_vs_matched_hcc_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_matched_hcc_tumor_vs_matched_hcc_normal.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#######
# set11

# 12/28/2020 note
# 8 samples, HWIDE recurrent
# 11 samples, HWIDE no-recurrent

# log2( mean(ix2) / mean (ix1) ) fold change
results_matched_hcc_HWIDE_recurrent_vs_hcc_HWIDE_no_recurrent_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=tumorSampleName_hcc_HWIDE_no_recurrent,ix2=tumorSampleName_hcc_HWIDE_recurrent)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

#x<-met
#ix1<-tumorSampleName_hcc_no_recurrent
#ix2<-tumorSampleName_hcc_HWIDE_recurrent

####
res<-results_matched_hcc_HWIDE_recurrent_vs_hcc_HWIDE_no_recurrent_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###

#res<-res[!(res$pathway %in% "Xenobiotics"),]
#res$padj2<-p.adjust(res$pvalue,method="BH")

#source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlotPublication_v3.R')

plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

#plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_hcc_HWIDE_recurrent_vs_hcc_HWIDE_no_recurrent_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_HWIDE_recurrent_vs_hcc_HWIDE_no_recurrent.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)


#######
# set0

# log2( mean(ix2) / mean (ix1) ) fold change
results_pd_ptctv_tumor_vs_pd_ptctv_normal_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=normalSampleName_pd_and_ptctv,ix2=tumorSampleName_pd_and_ptctv)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

####
res<-results_pd_ptctv_tumor_vs_pd_ptctv_normal_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###
plot<-makeVolcanoPlotPublication(res=res,
                                 pvalueCol = 'padj',
                                 threshold_pvalue = 0.05,
                                 threshold_log2fc = 1,
                                 numOfCandidateLabels = 5,
                                 show_legend_position = "none")

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_pd_ptctv_tumor_vs_pd_ptctv_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_pd_ptctv_tumor_vs_pd_ptctv_normal.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)


#######
# set0

normalSampleName_nonHCC<-c(normalSampleName_ha,normalSampleName_pd_and_ptctv)
tumorSampleName_nonHCC<-c(tumorSampleName_ha,tumorSampleName_pd_and_ptctv)

# log2( mean(ix2) / mean (ix1) ) fold change
results_ha_pd_ptctv_tumor_vs_ha_pd_ptctv_normal_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=normalSampleName_nonHCC,ix2=tumorSampleName_nonHCC)
#results_tumor_vs_normal_wilcox_test$metaboliteName<-rownames(results_tumor_vs_normal_wilcox_test)

####
res<-results_ha_pd_ptctv_tumor_vs_ha_pd_ptctv_normal_wilcox_test
res$name<-rownames(res)
res$pathway<-metinfo[rownames(res),'SUPER_PATHWAY']
###
plot<-makeVolcanoPlotPublication(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
#dir.create(filePath,recursive = TRUE)
fileName<-"HCC_Elisa_PQN_normalized_ha_pd_ptctv_tumor_vs_ha_pd_ptctv_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5 ,width=3.5)
print(plot)
dev.off()

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_ha_pd_ptctv_tumor_vs_ha_pd_ptctv_normal.txt"
fileName<-file.path(filePath,fileName)

write.table(res,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)
