library(cola)

######
# load name mapping from histogy sample name to metabolomics sample name
filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/hcc_cancercell"
fileName<-"clin_to_YD_sampleName.csv"
fileName<-file.path(filePath,fileName)

nameMapping<-read.table(file=fileName,header=TRUE,sep=",",stringsAsFactors = FALSE)
nameMapping<-nameMapping[(nchar(nameMapping$name)!=0),]
colnames(nameMapping)<-c("Sample","Sample_HCC")

nameMapping<-nameMapping[grepl("T$",nameMapping$Sample),]

######

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

clingen_subset<-clingen[,c("RTK.RAS.RAF.OR.PIK3.AKT.MTOR.alteration",
                           "TERT.DAXX.mutation",
                           "WCD.Chr7",
                           "LOH..haploidy.or.UPD.",
                           "Complex..mutation",
                           "Nuclear.mt.mutation",
                           "mtDNA")]
colnames(clingen_subset)<-c("MTORpathway","TERTmut","WCDChr7","LOHorUPD","mtDNA_complex","Nuclear_mtDNA","mtDNA")

clingen_subset$sampleName<-rownames(clingen_subset)

clingen_subset[!is.na(clingen_subset$mtDNA_complex),]$mtDNA_complex<-"Y"
clingen_subset[is.na(clingen_subset$mtDNA_complex),]$mtDNA_complex<-"N"

clingen_subset[!is.na(clingen_subset$Nuclear_mtDNA),]$Nuclear_mtDNA<-"Y"
clingen_subset[is.na(clingen_subset$Nuclear_mtDNA),]$Nuclear_mtDNA<-"N"

clingen_subset[!(clingen_subset$mtDNA=="WT"),]$mtDNA<-"Y"
clingen_subset[clingen_subset$mtDNA=="WT",]$mtDNA<-"N"


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

sampleInfo_subset<-sampleInfo[sampleInfo$CLIENT_IDENTIFIER %in% nameMapping$Sample_HCC,]
sampleInfo_subset<-sampleInfo_subset[nameMapping$Sample_HCC,]

######


filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/02_17_2021/metabolomics_HCC_tumor_32_samples/hclust_custom_with_YD9_outlier_v3"
fileName<-"metabolomics_HCC_tumor_32_samples.Rd"
fileName<-file.path(filePath,fileName)

load(file=fileName)

metabolomics_clusters_consensus<-metabolite_cc_clusters

res<-metabolomics_clusters_consensus["MAD:hclust_custom"]

set.seed(123)
get_signatures(res, k = 3)

#lt = functional_enrichment(res, k = 3)

dimension_reduction(res, method="PCA",k=3)
dimension_reduction(res, method="UMAP",k=3)

collect_classes(res,show_row_names=TRUE)

met_cluster_label<-cbind(get_classes(res, k = 3), get_membership(res, k = 3))


met_cluster_label_common<-met_cluster_label[rownames(met_cluster_label) %in% nameMapping$Sample_HCC,]
met_cluster_label_common<-met_cluster_label_common[nameMapping$Sample_HCC,]

met_matrixOfcluster<-model.matrix(~as.factor(met_cluster_label_common$class)-1)
rownames(met_matrixOfcluster)<-rownames(met_cluster_label_common)
colnames(met_matrixOfcluster)<-c("m1","m2","m3")


#####

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/02_17_2021/mRNA_HCC_tumor_53_samples/hclust_custom"
fileName<-"mRNA_HCC_tumor_53_samples.Rd"
fileName<-file.path(filePath,fileName)

load(file=fileName)

mRNA_clusters_consensus<-mRNA_cc_clusters

res<-mRNA_clusters_consensus["MAD:hclust_custom"]

set.seed(123)
get_signatures(res, k = 3)

#lt = functional_enrichment(res, k = 3)

dimension_reduction(res, method="PCA",k=3)
dimension_reduction(res, method="UMAP",k=3)

collect_classes(res,show_row_names=TRUE)

mRNA_cluster_label<-cbind(get_classes(res, k = 3), get_membership(res, k = 3))

mRNA_cluster_label_common<-mRNA_cluster_label[rownames(mRNA_cluster_label) %in% nameMapping$Sample,]
mRNA_cluster_label_common<-mRNA_cluster_label_common[nameMapping$Sample,]
rownames(mRNA_cluster_label_common)<-nameMapping$Sample_HCC

mRNA_matrixOfcluster<-model.matrix(~as.factor(mRNA_cluster_label_common$class)-1)
rownames(mRNA_matrixOfcluster)<-rownames(mRNA_cluster_label_common)
colnames(mRNA_matrixOfcluster)<-c("g1","g2","g3")


#######





#######

matrixOfclusters<-cbind(met_matrixOfcluster,mRNA_matrixOfcluster)

#matrixOfclusters<-cbind(met_matrixOfcluster,mRNA_matrixOfcluster,decon_matrixOfcluster)

#####


######

my_histology_col<-data.frame(histology=sampleInfo_subset$HISTOLOGY2,stringsAsFactors = FALSE)
rownames(my_histology_col)<-sampleInfo_subset$CLIENT_IDENTIFIER

my_genetic_col<-data.frame(clingen_subset[,c(1:7)],stringsAsFactors = FALSE)
my_genetic_col<-my_genetic_col[sampleInfo_subset$CLIENT_IDENTIFIER,]

my_met_cluster_col<-data.frame(met_cluster=met_cluster_label_common$class,stringsAsFactors = FALSE)
rownames(my_met_cluster_col)<-rownames(met_cluster_label_common)

my_mRNA_cluster_col<-data.frame(mRNA_cluster=mRNA_cluster_label_common$class,stringsAsFactors = FALSE)
rownames(my_mRNA_cluster_col)<-rownames(mRNA_cluster_label_common)

my_annotation<-cbind(my_histology_col,my_genetic_col,my_met_cluster_col,my_mRNA_cluster_col)




register_partition_methods(
  hclust_custom = function(mat, k, ...) {
    d<-dist(t(mat),method="euclidean")
    hc_result<-hclust(d,method="ward.D2")
    output<-cutree(hc_result, k = k)
    return(output)
  }
)

register_top_value_methods(
  originalOrder = function(mat,...) {
    
    constantVector<-rep(1,nrow(mat))
    names(constantVector)<-rownames(mat)
    return(constantVector)
  }
  
) 


my_annotation_color<-list(histology=c("HWIDE"="red","HMIN"="royalblue","HWIDE Recurrent"="green"),
                          MTORpathway=c("Y"="yellow","N"="grey50"),
                          TERTmut=c("Y"="yellow","N"="grey50"),
                          WCDChr7=c("Y"="yellow","N"="grey50","NA"="black"),
                          LOHorUPD=c("Y"="yellow","N"="grey50","NA"="black"),
                          mtDNA_complex=c("Y"="yellow","N"="grey50"),
                          Nuclear_mtDNA=c("Y"="yellow","N"="grey50"),
                          mtDNA=c("Y"="yellow","N"="grey50"),
                          met_cluster=c("1"="red","2"="royalblue","3"="yellow"),
                          mRNA_cluster=c("1"="red","2"="royalblue","3"="yellow"))

if(TRUE){

set.seed(123)  
res_ccp = consensus_partition(t(matrixOfclusters),
                          top_value_method = "originalOrder",
                          top_n = c(4,5),
                          #top_n = c(4,6,8),
                          partition_method = "hclust_custom",
                          max_k = 6,
                          sample_by = "row",
                          p_sampling = 0.8,
                          partition_repeat = 50,
                          scale_rows=FALSE,
                          anno = my_annotation,
                          anno_col = my_annotation_color,
                          verbose=TRUE)

select_partition_number(res_ccp)
get_stats(res_ccp)
suggest_best_k(res_ccp)

collect_classes(res_ccp,show_row_names=TRUE)


met_mRNA_cluster_label<-cbind(get_classes(res_ccp, k = 4), get_membership(res_ccp, k = 4))
consensus_sample_order<-rownames(met_mRNA_cluster_label[order(met_mRNA_cluster_label$class,
                                                     met_mRNA_cluster_label$silhouette),])

}


#outputDir<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/02_17_2021/combine_mRNA_met_cluster"
#dir.create(outputDir,recursive = TRUE)
#cola_report(res_ccp, output_dir = outputDir)

if(FALSE){
get_stats(res_ccp)
suggest_best_k(res_ccp)

collect_classes(res_ccp,show_row_names=TRUE)

combined_cluster_label<-cbind(get_classes(res_ccp, k = 3), get_membership(res_ccp, k = 3))

combined_cluster_label_common<-combined_cluster_label[nameMapping$Sample_HCC,]



check_label<-data.frame(met=met_cluster_label_common[,1],
                        mRNA=mRNA_cluster_label_common[,1],
                        consensus=combined_cluster_label_common[,1],
                        stringsAsFactors = FALSE)
rownames(check_label)<-nameMapping$Sample_HCC
}

#####
# check TME status
#####

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/hcc_cancercell"
#fileName<-"MergedDeconvolution.HCC.csv"
fileName<-"MergedDeconvolution_add_charoentong_set_02_15_2021.HCC.csv"
fileName<-file.path(filePath,fileName)  

# Grab the RNA deconvolution data
decon = read.csv(fileName,header = TRUE,row.names = 1)

decon_subset<-decon[nameMapping$Sample,]
rownames(decon_subset)<-nameMapping$Sample_HCC

#####

newMet<-t(decon_subset)

rowSD<-apply(newMet,1,sd)

#sdCutoff<-quantile(sort(rowSD,decreasing = TRUE),percentileThreshold)
sdCutoff<-0
rowSubset<-names(rowSD[rowSD>sdCutoff])

numberOfRowsLeft<-length(rowSubset)
numberOfRowsTotal<-length(rownames(newMet))

message(sprintf("Use %s / %s metabolites in the data",numberOfRowsLeft,numberOfRowsTotal))

newMetPCA<-newMet[rowSubset,]
dim(newMetPCA)

decon_subset<-t(newMetPCA)

#####

decon_immune_subset<-decon_subset[,c(1:95)]
decon_hallmark_subset<-decon_subset[,c(96:145)]
decon_kegg_subset<-decon_subset[,c(146:210)]


####
# check xCell score
####

if(FALSE){
  
filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/hcc_cancercell"
#fileName<-"MergedDeconvolution.HCC.csv"
fileName<-"xCell_02_15_2021.HCC.csv"
fileName<-file.path(filePath,fileName)  

# Grab the RNA deconvolution data
decon = read.csv(fileName,header = TRUE,row.names = 1)

decon<-t(decon)

decon_subset<-decon[nameMapping$Sample,]
rownames(decon_subset)<-nameMapping$Sample_HCC

}

#####

#decon$log2CYT<-log2(decon$CYT)

#log2CYT<-decon$log2CYT

#decon<-cbind(log2CYT,decon[,c(2:7)])
#decon<-cbind(log2CYT,decon[,c(2:7)])
#decon<-cbind(log2CYT,decon[,c(2:4)])
#decon<-cbind(log2CYT,decon[,c(2:43,45:94)])

if(FALSE){
  
  # 36 TME signatures
  selected_TME_signaures<-c("log2CYT",
                            "StromalScore","ImmuneScore","ESTIMATEScore",
                            "IIS","TIS","APM1",
                            "CD8.T.cells","T.cells","T.helper.cells",
                            "Tcm.cells","Tem.cells","Th1.cells",
                            "Th17.cells","Th2.cells","Treg.cells",
                            "aDC","B.cells","Cytotoxic.cells",
                            "DC","Eosinophils","iDC",
                            "Macrophages","Mast.cells","Neutrophils",
                            "NK.CD56bright.cells","NK.CD56dim.cells","NK.cells",
                            "pDC","Tfh.cells","Tgd.cells",
                            "PD1","PDL1","CTLA4","APM2","Angiogenesis"
  )
  
}

if(FALSE){

# 13 Bindea et al + 3 Charoentong et al TME signatures
selected_TME_signaures<-c(#"CYT",
  #"StromalScore","ImmuneScore","ESTIMATEScore",
  #"IIS","TIS","APM1",
  #"CD8.T.cells",
  #"T.cells",
  "T.helper.cells",
  "Tcm.cells","Tem.cells",
  #"Th1.cells","Th17.cells","Th2.cells",
  #"Treg.cells",
  "aDC","B.cells",
  #"Cytotoxic.cells",
  #"DC",
  "Eosinophils","iDC",
  "Macrophages","Mast.cells","Neutrophils",
  "NK.CD56bright.cells","NK.CD56dim.cells",
  #"NK.cells",
  #"pDC",
  "Tfh.cells",
  #"Tgd.cells",
  #"PD1","PDL1","CTLA4","APM2","Angiogenesis"
  "Tgd","Treg","CD8Tcells"
)

}

if(FALSE){
  selected_TME_signaures<-c(#"CYT",
    "log2CYT",
    "StromalScore","ImmuneScore","ESTIMATEScore",
    "IIS","TIS","APM1",
    "CD8.T.cells","T.cells","T.helper.cells",
    "Tcm.cells","Tem.cells","Th1.cells",
    "Th17.cells","Th2.cells","Treg.cells",
    "aDC","B.cells","Cytotoxic.cells",
    "DC","Eosinophils","iDC",
    "Macrophages","Mast.cells","Neutrophils",
    "NK.CD56bright.cells","NK.CD56dim.cells","NK.cells",
    "pDC","Tfh.cells","Tgd.cells",
    "PD1","PDL1","CTLA4","APM2","Angiogenesis",
    "Tgd","Treg","CD8Tcells"
  )
  
}



#decon_TME_subset<-decon[,selected_TME_signaures]

#decon<-decon_TME_subset

#decon<-t(decon_TME_subset)

######

decon_immune_subset<-t(decon_immune_subset)
decon_hallmark_subset<-t(decon_hallmark_subset)
decon_kegg_subset<-t(decon_kegg_subset)

######

####
# metabolite feature selection section
##

source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/diffMetaboliteAbundance.R')
source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlot_v3.R')

met_mRNA_cluster_label_ordered<-met_mRNA_cluster_label[consensus_sample_order,]

sample_class_1<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==1,])
sample_class_2<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==2,])
sample_class_3<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==3,])
sample_class_4<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==4,])

#####

results_matched_hcc_tumor_class_1_vs_other_wilcox_test<-diffMetaboliteAbundance(x=decon_immune_subset,ix1=c(sample_class_2,sample_class_3,sample_class_4),ix2=sample_class_1)
res_wilcox_class_1<-results_matched_hcc_tumor_class_1_vs_other_wilcox_test
res_wilcox_class_1$name<-rownames(res_wilcox_class_1)
#res_wilcox_class_1$pathway<-metinfo[rownames(res_wilcox_class_1),'SUPER_PATHWAY']
plot<-makeVolcanoPlot(res=res_wilcox_class_1,pvalueCol = 'padj',threshold_pvalue = 0.1,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/06_28_2021/combine_mRNA_met_cluster/TME"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_integrated_clustering_tumor_class_1_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/06_28_2021/combine_mRNA_met_cluster/TME"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_class_1_vs_others.txt"
fileName<-file.path(filePath,fileName)

write.table(res_wilcox_class_1,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)


#####

results_matched_hcc_tumor_class_2_vs_other_wilcox_test<-diffMetaboliteAbundance(x=decon_immune_subset,ix1=c(sample_class_1,sample_class_3,sample_class_4),ix2=sample_class_2)
res_wilcox_class_2<-results_matched_hcc_tumor_class_2_vs_other_wilcox_test
res_wilcox_class_2$name<-rownames(res_wilcox_class_2)
#res_wilcox_class_2$pathway<-metinfo[rownames(res_wilcox_class_2),'SUPER_PATHWAY']
plot<-makeVolcanoPlot(res=res_wilcox_class_2,pvalueCol = 'padj',threshold_pvalue = 0.1,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/06_28_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_integrated_clustering_tumor_class_2_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/06_28_2021/combine_mRNA_met_cluster/metabolomics"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_class_2_vs_others.txt"
fileName<-file.path(filePath,fileName)

write.table(res_wilcox_class_2,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#####

results_matched_hcc_tumor_class_3_vs_other_wilcox_test<-diffMetaboliteAbundance(x=decon_immune_subset,ix1=c(sample_class_1,sample_class_2,sample_class_4),ix2=sample_class_3)
res_wilcox_class_3<-results_matched_hcc_tumor_class_3_vs_other_wilcox_test
res_wilcox_class_3$name<-rownames(res_wilcox_class_3)
#res_wilcox_class_3$pathway<-metinfo[rownames(res_wilcox_class_3),'SUPER_PATHWAY']
plot<-makeVolcanoPlot(res=res_wilcox_class_3,pvalueCol = 'padj',threshold_pvalue = 0.1,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/06_28_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_integrated_clustering_tumor_class_3_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/06_28_2021/combine_mRNA_met_cluster/metabolomics"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_class_3_vs_others.txt"
fileName<-file.path(filePath,fileName)

write.table(res_wilcox_class_3,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#####

results_matched_hcc_tumor_class_4_vs_other_wilcox_test<-diffMetaboliteAbundance(x=decon_immune_subset,ix1=c(sample_class_1,sample_class_2,sample_class_3),ix2=sample_class_4)
res_wilcox_class_4<-results_matched_hcc_tumor_class_4_vs_other_wilcox_test
res_wilcox_class_4$name<-rownames(res_wilcox_class_4)
#res_wilcox_class_4$pathway<-metinfo[rownames(res_wilcox_class_4),'SUPER_PATHWAY']
plot<-makeVolcanoPlot(res=res_wilcox_class_1,pvalueCol = 'padj',threshold_pvalue = 0.1,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/06_28_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_integrated_clustering_tumor_class_4_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/06_28_2021/combine_mRNA_met_cluster/metabolomics"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_class_4_vs_others.txt"
fileName<-file.path(filePath,fileName)

write.table(res_wilcox_class_4,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#####

set1<-rownames(results_matched_hcc_tumor_class_1_vs_other_wilcox_test[results_matched_hcc_tumor_class_1_vs_other_wilcox_test$padj<0.1,])
set2<-rownames(results_matched_hcc_tumor_class_2_vs_other_wilcox_test[results_matched_hcc_tumor_class_2_vs_other_wilcox_test$padj<0.1,])
set3<-rownames(results_matched_hcc_tumor_class_3_vs_other_wilcox_test[results_matched_hcc_tumor_class_3_vs_other_wilcox_test$padj<0.1,])
set4<-rownames(results_matched_hcc_tumor_class_3_vs_other_wilcox_test[results_matched_hcc_tumor_class_3_vs_other_wilcox_test$padj<0.1,])

set1<-rownames(results_matched_hcc_tumor_class_1_vs_other_wilcox_test[results_matched_hcc_tumor_class_1_vs_other_wilcox_test$pvalue<0.1,])
set2<-rownames(results_matched_hcc_tumor_class_2_vs_other_wilcox_test[results_matched_hcc_tumor_class_2_vs_other_wilcox_test$pvalue<0.1,])
set3<-rownames(results_matched_hcc_tumor_class_3_vs_other_wilcox_test[results_matched_hcc_tumor_class_3_vs_other_wilcox_test$pvalue<0.1,])
set4<-rownames(results_matched_hcc_tumor_class_3_vs_other_wilcox_test[results_matched_hcc_tumor_class_3_vs_other_wilcox_test$pvalue<0.1,])


significant_tme_set<-unique(c(set1,
                                     set2,
                                     set3,
                                     set4))


message(sprintf("%s sigfinicant candidates",length(significant_tme_set)))


significant_tme_set2<-significant_tme_set[significant_tme_set %in% selected_TME_signaures]

message(sprintf("%s sigfinicant candidates",length(significant_tme_set2)))

######
library(pheatmap)

#decon<-get_matrix(res)

#decon_common<-decon_immune_subset[,nameMapping$Sample]
#colnames(decon_common)<-nameMapping$Sample_HCC

decon_common<-decon_immune_subset

decon_common<-decon_common[,consensus_sample_order]

my_consensus_cluster_col<-data.frame(consensus_cluster=met_mRNA_cluster_label$class,stringsAsFactors = FALSE)
rownames(my_consensus_cluster_col)<-rownames(met_mRNA_cluster_label)

my_annotation2<-cbind(my_consensus_cluster_col,my_annotation)
my_annotation2[is.na(my_annotation2)]<-"N/A"
my_annotation2<-my_annotation2[consensus_sample_order,]

my_annotation2$mRNA_cluster<-as.character(my_annotation2$mRNA_cluster)
my_annotation2$met_cluster<-as.character(my_annotation2$met_cluster)
my_annotation2$consensus_cluster<-as.character(my_annotation2$consensus_cluster)

my_annotation_color2<-list(histology=c("HWIDE"="indianred1","HMIN"="skyblue1","HWIDE Recurrent"="green"),
                           MTORpathway=c("Y"="yellow","N"="grey50"),
                           TERTmut=c("Y"="yellow","N"="grey50"),
                           WCDChr7=c("Y"="yellow","N"="grey50","N/A"="grey90"),
                           LOHorUPD=c("Y"="yellow","N"="grey50","N/A"="grey90"),
                           mtDNA_complex=c("Y"="yellow","N"="grey50"),
                           Nuclear_mtDNA=c("Y"="yellow","N"="grey50"),
                           mtDNA=c("Y"="yellow","N"="grey50"),
                           met_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1"),
                           mRNA_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1"),
                           #decon_cluster=c("1"="red","2"="royalblue","3"="yellow"),
                           #decon_cluster=c("1"="red","2"="royalblue","3"="yellow","4"="green","5"="pink"),
                           consensus_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1","4"="pink1"))

data_subset<-decon_common

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(data_subset, 1, cal_z_score))

#data_subset_norm[data_subset_norm>3]<-3
#data_subset_norm[data_subset_norm<(-3)]<-(-3)

data_subset_norm[data_subset_norm>2]<-2
data_subset_norm[data_subset_norm<(-2)]<-(-2)

#data_subset_norm[data_subset_norm>1]<-1
#data_subset_norm[data_subset_norm<(-1)]<-(-1)

paletteLength <- 15
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(data_subset_norm, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(data_subset_norm, na.rm=TRUE)/paletteLength, max(data_subset_norm, na.rm=TRUE), length.out=floor(paletteLength/2)))

graph<-pheatmap(data_subset_norm[significant_tme_set2,],
         color =myColor, 
         breaks=myBreaks,
         gaps_col = c(7,16,22),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_col = my_annotation2,
         annotation_colors = my_annotation_color2,
         scale="none",
         fontsize = 8,
         border_color = "black")

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/02_19_2021/combine_mRNA_met_cluster/TME"
dir.create(filePath,recursive = TRUE)
fileName<-paste("HCC_tme_cc_cluster_full_immune_signatures_set.pdf",sep="")
fileName<-file.path(filePath,fileName)

pdf(fileName,width=9,height=15)
print(graph)
dev.off()

#####
# check on hallmark set
#####

decon_common<-decon_hallmark_subset

decon_common<-decon_common[,consensus_sample_order]

my_consensus_cluster_col<-data.frame(consensus_cluster=met_mRNA_cluster_label$class,stringsAsFactors = FALSE)
rownames(my_consensus_cluster_col)<-rownames(met_mRNA_cluster_label)

my_annotation2<-cbind(my_consensus_cluster_col,my_annotation)
my_annotation2[is.na(my_annotation2)]<-"N/A"
my_annotation2<-my_annotation2[consensus_sample_order,]

my_annotation2$mRNA_cluster<-as.character(my_annotation2$mRNA_cluster)
my_annotation2$met_cluster<-as.character(my_annotation2$met_cluster)
my_annotation2$consensus_cluster<-as.character(my_annotation2$consensus_cluster)

my_annotation_color2<-list(histology=c("HWIDE"="indianred1","HMIN"="skyblue1","HWIDE Recurrent"="green"),
                           MTORpathway=c("Y"="yellow","N"="grey50"),
                           TERTmut=c("Y"="yellow","N"="grey50"),
                           WCDChr7=c("Y"="yellow","N"="grey50","N/A"="grey90"),
                           LOHorUPD=c("Y"="yellow","N"="grey50","N/A"="grey90"),
                           mtDNA_complex=c("Y"="yellow","N"="grey50"),
                           Nuclear_mtDNA=c("Y"="yellow","N"="grey50"),
                           mtDNA=c("Y"="yellow","N"="grey50"),
                           met_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1"),
                           mRNA_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1"),
                           #decon_cluster=c("1"="red","2"="royalblue","3"="yellow"),
                           #decon_cluster=c("1"="red","2"="royalblue","3"="yellow","4"="green","5"="pink"),
                           consensus_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1","4"="pink1"))

data_subset<-decon_common

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(data_subset, 1, cal_z_score))

#data_subset_norm[data_subset_norm>3]<-3
#data_subset_norm[data_subset_norm<(-3)]<-(-3)

data_subset_norm[data_subset_norm>2]<-2
data_subset_norm[data_subset_norm<(-2)]<-(-2)

#data_subset_norm[data_subset_norm>1]<-1
#data_subset_norm[data_subset_norm<(-1)]<-(-1)

paletteLength <- 15
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(data_subset_norm, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(data_subset_norm, na.rm=TRUE)/paletteLength, max(data_subset_norm, na.rm=TRUE), length.out=floor(paletteLength/2)))

graph<-pheatmap(data_subset_norm,
                color =myColor, 
                breaks=myBreaks,
                gaps_col = c(7,16,22),
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                annotation_col = my_annotation2,
                annotation_colors = my_annotation_color2,
                scale="none",
                fontsize = 8,
                border_color = "black")

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/02_19_2021/combine_mRNA_met_cluster/TME"
dir.create(filePath,recursive = TRUE)
fileName<-paste("HCC_tme_cc_cluster_hallmark_signatures_set.pdf",sep="")
fileName<-file.path(filePath,fileName)

pdf(fileName,width=12,height=15)
print(graph)
dev.off()

#####
# check KEGG pathways
#####

decon_common<-decon_kegg_subset

decon_common<-decon_common[,consensus_sample_order]

my_consensus_cluster_col<-data.frame(consensus_cluster=met_mRNA_cluster_label$class,stringsAsFactors = FALSE)
rownames(my_consensus_cluster_col)<-rownames(met_mRNA_cluster_label)

my_annotation2<-cbind(my_consensus_cluster_col,my_annotation)
my_annotation2[is.na(my_annotation2)]<-"N/A"
my_annotation2<-my_annotation2[consensus_sample_order,]

my_annotation2$mRNA_cluster<-as.character(my_annotation2$mRNA_cluster)
my_annotation2$met_cluster<-as.character(my_annotation2$met_cluster)
my_annotation2$consensus_cluster<-as.character(my_annotation2$consensus_cluster)

my_annotation_color2<-list(histology=c("HWIDE"="indianred1","HMIN"="skyblue1","HWIDE Recurrent"="green"),
                           MTORpathway=c("Y"="yellow","N"="grey50"),
                           TERTmut=c("Y"="yellow","N"="grey50"),
                           WCDChr7=c("Y"="yellow","N"="grey50","N/A"="grey90"),
                           LOHorUPD=c("Y"="yellow","N"="grey50","N/A"="grey90"),
                           mtDNA_complex=c("Y"="yellow","N"="grey50"),
                           Nuclear_mtDNA=c("Y"="yellow","N"="grey50"),
                           mtDNA=c("Y"="yellow","N"="grey50"),
                           met_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1"),
                           mRNA_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1"),
                           #decon_cluster=c("1"="red","2"="royalblue","3"="yellow"),
                           #decon_cluster=c("1"="red","2"="royalblue","3"="yellow","4"="green","5"="pink"),
                           consensus_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1","4"="pink1"))

data_subset<-decon_common

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(data_subset, 1, cal_z_score))

#data_subset_norm[data_subset_norm>3]<-3
#data_subset_norm[data_subset_norm<(-3)]<-(-3)

data_subset_norm[data_subset_norm>2]<-2
data_subset_norm[data_subset_norm<(-2)]<-(-2)

#data_subset_norm[data_subset_norm>1]<-1
#data_subset_norm[data_subset_norm<(-1)]<-(-1)

paletteLength <- 15
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(data_subset_norm, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(data_subset_norm, na.rm=TRUE)/paletteLength, max(data_subset_norm, na.rm=TRUE), length.out=floor(paletteLength/2)))

graph<-pheatmap(data_subset_norm,
                color =myColor, 
                breaks=myBreaks,
                gaps_col = c(7,16,22),
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                annotation_col = my_annotation2,
                annotation_colors = my_annotation_color2,
                scale="none",
                fontsize = 8,
                border_color = "black")

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/02_19_2021/combine_mRNA_met_cluster/TME"
dir.create(filePath,recursive = TRUE)
fileName<-paste("HCC_tme_cc_cluster_KEGG_signatures_set.pdf",sep="")
fileName<-file.path(filePath,fileName)

pdf(fileName,width=12,height=15)
print(graph)
dev.off()
