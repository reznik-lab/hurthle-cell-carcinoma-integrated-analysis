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

#######
# additional annotation correction 03/25/2021
######


clingen["YD14","RTK.RAS.RAF.OR.PIK3.AKT.MTOR.alteration"]<-"Y"
clingen["YD14","WCD.Chr7"]<-"Y"

clingen["YD41","RTK.RAS.RAF.OR.PIK3.AKT.MTOR.alteration"]<-"N"


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

metinfoSubset<-metinfo[,c("name","H_KEGG")]

#####


newMet<-met[,nameMapping$Sample_HCC]

rowSD<-apply(newMet,1,sd)

#sdCutoff<-quantile(sort(rowSD,decreasing = TRUE),percentileThreshold)
sdCutoff<-0
rowSubset<-names(rowSD[rowSD>sdCutoff])

numberOfRowsLeft<-length(rowSubset)
numberOfRowsTotal<-length(rownames(newMet))

message(sprintf("Use %s / %s metabolites in the data",numberOfRowsLeft,numberOfRowsTotal))

newMetPCA<-newMet[rowSubset,]
dim(newMetPCA)

met_tumor<-newMetPCA

log2_met_tumor<-log2(met_tumor)




######
# load coca results
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
                          #top_n = c(3,4,5),
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

sampleName<-rownames(get_classes(res_ccp, k = 4))

met_mRNA_cluster_label<-cbind(sampleName,get_classes(res_ccp, k = 4), get_membership(res_ccp, k = 4),get_anno(res_ccp))
consensus_sample_order<-rownames(met_mRNA_cluster_label[order(met_mRNA_cluster_label$class,
                                                     met_mRNA_cluster_label$silhouette),])

met_mRNA_cluster_label<-met_mRNA_cluster_label[consensus_sample_order,]

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster"
dir.create(filePath,recursive = TRUE)
fileName<-"met_mRNA_consensus_cluster_assignment.txt"
fileName<-file.path(filePath,fileName)

write.table(met_mRNA_cluster_label,file=fileName,sep="\t",quote=FALSE,row.names = FALSE, col.names = TRUE)


}

outputDir<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster"
dir.create(outputDir,recursive = TRUE)
fileName<-"mRNA_met_coca_results.Rd"
fileName<-file.path(outputDir,fileName)

save(res_ccp,file=fileName)


outputDir<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster"
dir.create(outputDir,recursive = TRUE)
cola_report(res_ccp, output_dir = outputDir)


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
# make metaboloite heatmap based on COCA
#####

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

results_matched_hcc_tumor_class_1_vs_other_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=c(sample_class_2,sample_class_3,sample_class_4),ix2=sample_class_1)
res_wilcox_class_1<-results_matched_hcc_tumor_class_1_vs_other_wilcox_test
res_wilcox_class_1$name<-rownames(res_wilcox_class_1)
res_wilcox_class_1$pathway<-metinfo[rownames(res_wilcox_class_1),'SUPER_PATHWAY']
plot<-makeVolcanoPlot(res=res_wilcox_class_1,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_integrated_clustering_tumor_class_1_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_class_1_vs_others.txt"
fileName<-file.path(filePath,fileName)

write.table(res_wilcox_class_1,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)


#####

results_matched_hcc_tumor_class_2_vs_other_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=c(sample_class_1,sample_class_3,sample_class_4),ix2=sample_class_2)
res_wilcox_class_2<-results_matched_hcc_tumor_class_2_vs_other_wilcox_test
res_wilcox_class_2$name<-rownames(res_wilcox_class_2)
res_wilcox_class_2$pathway<-metinfo[rownames(res_wilcox_class_2),'SUPER_PATHWAY']
plot<-makeVolcanoPlot(res=res_wilcox_class_2,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_integrated_clustering_tumor_class_2_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_class_2_vs_others.txt"
fileName<-file.path(filePath,fileName)

write.table(res_wilcox_class_2,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#####

results_matched_hcc_tumor_class_3_vs_other_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=c(sample_class_1,sample_class_2,sample_class_4),ix2=sample_class_3)
res_wilcox_class_3<-results_matched_hcc_tumor_class_3_vs_other_wilcox_test
res_wilcox_class_3$name<-rownames(res_wilcox_class_3)
res_wilcox_class_3$pathway<-metinfo[rownames(res_wilcox_class_3),'SUPER_PATHWAY']
plot<-makeVolcanoPlot(res=res_wilcox_class_3,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_integrated_clustering_tumor_class_3_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_class_3_vs_others.txt"
fileName<-file.path(filePath,fileName)

write.table(res_wilcox_class_3,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#####

results_matched_hcc_tumor_class_4_vs_other_wilcox_test<-diffMetaboliteAbundance(x=met,ix1=c(sample_class_1,sample_class_2,sample_class_3),ix2=sample_class_4)
res_wilcox_class_4<-results_matched_hcc_tumor_class_4_vs_other_wilcox_test
res_wilcox_class_4$name<-rownames(res_wilcox_class_4)
res_wilcox_class_4$pathway<-metinfo[rownames(res_wilcox_class_4),'SUPER_PATHWAY']
plot<-makeVolcanoPlot(res=res_wilcox_class_4,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 1)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_integrated_clustering_tumor_class_4_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_class_4_vs_others.txt"
fileName<-file.path(filePath,fileName)

write.table(res_wilcox_class_4,file=fileName,sep="\t",quote=FALSE,row.names=FALSE,col.names = TRUE)

#####

padj_threshold<-0.05
log2FC_threshold<-1

set1<-rownames(res_wilcox_class_1[res_wilcox_class_1$padj<padj_threshold & abs(res_wilcox_class_1$log2FC)> log2FC_threshold,])
set2<-rownames(res_wilcox_class_2[res_wilcox_class_2$padj<padj_threshold & abs(res_wilcox_class_2$log2FC)> log2FC_threshold,])
set3<-rownames(res_wilcox_class_3[res_wilcox_class_3$padj<padj_threshold & abs(res_wilcox_class_3$log2FC)> log2FC_threshold,])
set4<-rownames(res_wilcox_class_4[res_wilcox_class_4$padj<padj_threshold & abs(res_wilcox_class_4$log2FC)> log2FC_threshold,])


significant_metabolite_set<-unique(c(set1,
                                     set2,
                                     set3,
                                     set4))


message(sprintf("%s sigfinicant candidates",length(significant_metabolite_set)))


######
#library(pheatmap)
library(ComplexHeatmap)

#decon<-get_matrix(res)

met_common<-log2_met_tumor

met_common<-met_common[,consensus_sample_order]

my_consensus_cluster_col<-data.frame(consensus_cluster=met_mRNA_cluster_label$class,stringsAsFactors = FALSE)
rownames(my_consensus_cluster_col)<-rownames(met_mRNA_cluster_label)

my_annotation<-my_annotation[consensus_sample_order,]

my_annotation2<-cbind(my_consensus_cluster_col,my_annotation)
my_annotation2[is.na(my_annotation2)]<-"N/A"

my_annotation2<-my_annotation2[consensus_sample_order,]


my_annotation_color2<-list(histology=c("HWIDE"="indianred1","HMIN"="skyblue1","HWIDE Recurrent"="green"),
                          MTORpathway=c("Y"="yellow","N"="grey50"),
                          TERTmut=c("Y"="yellow","N"="grey50"),
                          WCDChr7=c("Y"="yellow","N"="grey50","N/A"="grey90"),
                          LOHorUPD=c("Y"="yellow","N"="grey50","N/A"="grey90"),
                          mtDNA_complex=c("Y"="yellow","N"="grey50"),
                          Nuclear_mtDNA=c("Y"="yellow","N"="grey50"),
                          mtDNA=c("Y"="yellow","N"="grey50"),
                          met_cluster=c("1"="mediumpurple1","2"="red","3"="royalblue"),
                          mRNA_cluster=c("1"="mediumpurple1","2"="red","3"="royalblue"),
                          #decon_cluster=c("1"="red","2"="royalblue","3"="yellow"),
                          #decon_cluster=c("1"="red","2"="royalblue","3"="yellow","4"="green","5"="pink"),
                          consensus_cluster=c("1"="pink1","2"="red","3"="royalblue","4"="mediumpurple1"))

data_subset<-met_common

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
myBreaks <- c(seq(min(data_subset_norm), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(data_subset_norm)/paletteLength, max(data_subset_norm), length.out=floor(paletteLength/2)))

additional_metabolite_set<-c("saccharopine", "2-aminoadipate", "3-hydroxyglutarate","citrate","aconitate [cis or trans]","cysteine","lactate")
additional_metabolite_set2<-c("lactate","quinolinate",
                              "oleate/vaccenate (18:1)","stearate (18:0)","palmitate (16:0)","histamine",
                              "nicotinamide adenine dinucleotide (NAD+)",
                              "orotate")

graph<-ComplexHeatmap::pheatmap(data_subset_norm[unique(c(additional_metabolite_set,additional_metabolite_set2,significant_metabolite_set)),],
         color =myColor, 
         breaks=myBreaks,
         #gaps_col = c(7,13,19),
         gaps_col = c(7,16,22),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_col = my_annotation2,
         annotation_colors = my_annotation_color2,
         scale="none",
         fontsize = 7,
         border_color = "black")

graph<-ComplexHeatmap::pheatmap(data_subset_norm[unique(c(additional_metabolite_set,significant_metabolite_set)),],
                                color =myColor, 
                                breaks=myBreaks,
                                #gaps_col = c(7,13,19),
                                gaps_col = c(7,16,22),
                                cluster_cols = FALSE,
                                cluster_rows = TRUE,
                                annotation_col = my_annotation2,
                                annotation_colors = my_annotation_color2,
                                scale="none",
                                fontsize = 7,
                                border_color = "black")

column_split_order<-met_mRNA_cluster_label_ordered$class

column_ha = HeatmapAnnotation(
  mRNA_cluster = my_annotation2$mRNA_cluster,
  met_cluster = my_annotation2$met_cluster,
  mtDNA = my_annotation2$mtDNA,
  Nuclear_mtDNA = my_annotation2$Nuclear_mtDNA,
  LOHorUPD = my_annotation2$LOHorUPD,
  WCDChr7 = my_annotation2$WCDChr7,
  TERTmut = my_annotation2$TERTmut,
  MTORpathway = my_annotation2$MTORpathway,
  histology = my_annotation2$histology,
  consensus_cluster = my_annotation2$consensus_cluster,
  col=my_annotation_color2,
  gp = gpar(col="black"),
  border=FALSE,
  annotation_name_gp = gpar(fontsize=7),
  gap = unit(0.5, "mm"),
  #annotation_height = unit(c(0.5, 0.5,0.5, 0.5,0.5, 0.5,0.5, 0.5,0.5, 0.5), "cm")
  simple_anno_size = unit(2.5,"mm"),
  annotation_legend_param = list(title_gp = gpar(fontsize = 7), 
                                 labels_gp = gpar(fontsize = 7),
                                 border = "black",
                                 grid_height = unit(2.5, "mm"),
                                 grid_width = unit(2.5, "mm"))
)


ht<-Heatmap(data_subset_norm[unique(c(additional_metabolite_set,significant_metabolite_set)),],
        name = "mat",
        rect_gp = gpar(col = "black", lwd = 1),
        top_annotation = column_ha,
        col = myColor,
        cluster_columns = FALSE,
        column_names_rot = 90,
        column_split = column_split_order,
        column_title = NULL,
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 7))


draw(ht, heatmap_legend_side="bottom", annotation_legend_side="bottom",merge_legend=TRUE)


filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-paste("HCC_metabolite_cc_cluster_sig_padj_0.05.pdf",sep="")
fileName<-file.path(filePath,fileName)

pdf(fileName,width=5,height=7.2)
draw(ht, heatmap_legend_side="bottom", annotation_legend_side="bottom",merge_legend=TRUE)
dev.off()


#####



#####

graph<-pheatmap(data_subset_norm[unique(c(additional_metabolite_set2,significant_metabolite_set)),],
                color =myColor, 
                breaks=myBreaks,
                #gaps_col = c(7,13,19),
                gaps_col = c(7,16,22),
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                annotation_col = my_annotation2,
                annotation_colors = my_annotation_color2,
                scale="none",
                fontsize = 8,
                border_color = "black")


filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-paste("HCC_metabolite_cc_cluster_sig_additional_set_2_padj_0.05.pdf",sep="")
fileName<-file.path(filePath,fileName)

pdf(fileName,width=9,height=12)
print(graph)
dev.off()

######

#########
# run GSEA
#########
library(fgsea)

met_mRNA_cluster_label_ordered<-met_mRNA_cluster_label[consensus_sample_order,]

sample_class_1<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==1,])
sample_class_2<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==2,])
sample_class_3<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==3,])
sample_class_4<-rownames(met_mRNA_cluster_label_ordered[met_mRNA_cluster_label_ordered$class==4,])

met_rowMeans<-rowMeans(log2_met_tumor)
rowMean_threshold<-quantile(met_rowMeans,0.05)
selectedMet1<-names(met_rowMeans[met_rowMeans>rowMean_threshold])

row_non_empty_rate<-apply(log2_met_tumor,1,function(x){sum(x!=0)/length(x)})
selectedMet2<-names(row_non_empty_rate[row_non_empty_rate>0.7])

#selectedMets<-intersect(selectedMet1,selectedMet2)
selectedMets<-selectedMet1

#selectedMets<-significant_metabolite_set

log2_met_tumor_filtered<-log2_met_tumor[rownames(log2_met_tumor) %in% selectedMets,]

#dat<-log2_met_tumor_filtered
#ix1<-sample_class_1

signalToNoiseRatio<-function(dat,ix1){
  
  output<-list()  
  
  output<-sapply( 1:nrow(dat), function(idx) {
    
    vector1<-dat[idx,ix1] 
    vector2<-dat[idx,!(colnames(dat) %in% ix1)]    
    
    uA<-mean(vector1)     
    uB<-mean(vector2)
    
    sdA<-sd(vector1)
    sdB<-sd(vector2)
    
    signalToNoiseRatio<-(uA-uB)/(sdA+sdB)
    
    output[[idx]]<-signalToNoiseRatio
    
  })
  
  names(output)<-rownames(dat)
  
  return(output)
}

signalToNoiseRatio_class_1<-signalToNoiseRatio(log2_met_tumor_filtered,sample_class_1)
signalToNoiseRatio_class_2<-signalToNoiseRatio(log2_met_tumor_filtered,sample_class_2)
signalToNoiseRatio_class_3<-signalToNoiseRatio(log2_met_tumor_filtered,sample_class_3)
signalToNoiseRatio_class_4<-signalToNoiseRatio(log2_met_tumor_filtered,sample_class_4)

####
# change metabolite name to KEGG ID
####
# Todo fix aggregate name
####

#vectorInput<-signalToNoiseRatio_KEGG_class_1

fixAggregateKEGGvector<-function(vectorInput){
  
      broken_set<-vectorInput[grepl(",",names(vectorInput))]
      good_set<-vectorInput[!grepl(",",names(vectorInput))]  

      str_list<-strsplit(names(broken_set),",")
      
      fixed_set<-list()
      #fixed_vector<-NULL
      for( idx in 1:length(broken_set) ){
        
        fixed_vector<-rep(broken_set[idx],length(str_list[[idx]]))
        names(fixed_vector)<-str_list[[idx]]
        fixed_set[[idx]]<-fixed_vector
        
      }  
      
      fixed_set<-unlist(fixed_set)
      fixed_set<-c(fixed_set,good_set)
      
      return(fixed_set)
}




metNameMapping<-metinfoSubset[metinfoSubset$name %in% names(signalToNoiseRatio_class_1),]
metNameMapping_KEGG<-metNameMapping[!is.na(metNameMapping$H_KEGG),]
signalToNoiseRatio_KEGG_class_1<-signalToNoiseRatio_class_1[metNameMapping_KEGG$name]
names(signalToNoiseRatio_KEGG_class_1)<-metNameMapping_KEGG$H_KEGG

signalToNoiseRatio_KEGG_class_1<-signalToNoiseRatio_KEGG_class_1[unique(names(signalToNoiseRatio_KEGG_class_1))]
signalToNoiseRatio_KEGG_class_1<-fixAggregateKEGGvector(signalToNoiseRatio_KEGG_class_1)

metNameMapping<-metinfoSubset[metinfoSubset$name %in% names(signalToNoiseRatio_class_2),]
metNameMapping_KEGG<-metNameMapping[!is.na(metNameMapping$H_KEGG),]
signalToNoiseRatio_KEGG_class_2<-signalToNoiseRatio_class_2[metNameMapping_KEGG$name]
names(signalToNoiseRatio_KEGG_class_2)<-metNameMapping_KEGG$H_KEGG

signalToNoiseRatio_KEGG_class_2<-signalToNoiseRatio_KEGG_class_2[unique(names(signalToNoiseRatio_KEGG_class_2))]
signalToNoiseRatio_KEGG_class_2<-fixAggregateKEGGvector(signalToNoiseRatio_KEGG_class_2)

metNameMapping<-metinfoSubset[metinfoSubset$name %in% names(signalToNoiseRatio_class_3),]
metNameMapping_KEGG<-metNameMapping[!is.na(metNameMapping$H_KEGG),]
signalToNoiseRatio_KEGG_class_3<-signalToNoiseRatio_class_3[metNameMapping_KEGG$name]
names(signalToNoiseRatio_KEGG_class_3)<-metNameMapping_KEGG$H_KEGG

signalToNoiseRatio_KEGG_class_3<-signalToNoiseRatio_KEGG_class_3[unique(names(signalToNoiseRatio_KEGG_class_3))]
signalToNoiseRatio_KEGG_class_3<-fixAggregateKEGGvector(signalToNoiseRatio_KEGG_class_3)

metNameMapping<-metinfoSubset[metinfoSubset$name %in% names(signalToNoiseRatio_class_4),]
metNameMapping_KEGG<-metNameMapping[!is.na(metNameMapping$H_KEGG),]
signalToNoiseRatio_KEGG_class_4<-signalToNoiseRatio_class_4[metNameMapping_KEGG$name]
names(signalToNoiseRatio_KEGG_class_4)<-metNameMapping_KEGG$H_KEGG

signalToNoiseRatio_KEGG_class_4<-signalToNoiseRatio_KEGG_class_4[unique(names(signalToNoiseRatio_KEGG_class_4))]
signalToNoiseRatio_KEGG_class_4<-fixAggregateKEGGvector(signalToNoiseRatio_KEGG_class_4)

####

abs_signalToNoiseRatio_rank_KEGG_class_1<-sort(abs(signalToNoiseRatio_KEGG_class_1),decreasing = TRUE)
abs_signalToNoiseRatio_rank_KEGG_class_2<-sort(abs(signalToNoiseRatio_KEGG_class_2),decreasing = TRUE)
abs_signalToNoiseRatio_rank_KEGG_class_3<-sort(abs(signalToNoiseRatio_KEGG_class_3),decreasing = TRUE)
abs_signalToNoiseRatio_rank_KEGG_class_4<-sort(abs(signalToNoiseRatio_KEGG_class_4),decreasing = TRUE)

signalToNoiseRatio_rank_KEGG_class_1<-sort(signalToNoiseRatio_KEGG_class_1,decreasing = TRUE)
signalToNoiseRatio_rank_KEGG_class_2<-sort(signalToNoiseRatio_KEGG_class_2,decreasing = TRUE)
signalToNoiseRatio_rank_KEGG_class_3<-sort(signalToNoiseRatio_KEGG_class_3,decreasing = TRUE)
signalToNoiseRatio_rank_KEGG_class_4<-sort(signalToNoiseRatio_KEGG_class_4,decreasing = TRUE)


#####

#######
# run enricher on differential metabolite abundance
#######
if(FALSE){
filePath<-"~/work/Ed_lab/HCC_project/data/KEGG_metabolic_pathway/parsed_data"
fileName<-paste("KEGG_metabolic_pathway_metaboliteInfoAll.txt",sep="")
fileName<-file.path(filePath,fileName)

dat<-fread(file=fileName,header=TRUE,sep="\t",stringsAsFactors = FALSE)

#####
library(clusterProfiler)

set.seed(42)

pathwayToMetabolite<-dat[,c("pathwayName","keggMetaboliteID")]

pathwayToMetaboliteName<-dat[,c("pathwayName","pathwayName")]

hh<-data.frame(names(signalToNoiseRatio_class_2),signalToNoiseRatio_class_2,stringsAsFactors = FALSE)
colnames(hh)<-c("name","snRatio")

matchedDiffMetAbundance<-merge(hh,metinfo,by="name",all.x=TRUE)
matchedMetAbundanceList<-matchedDiffMetAbundance$H_KEGG

x<-enricher(matchedMetAbundanceList,
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            minGSSize = 3,
            maxGSSize = 500,
            qvalueCutoff = 1,
            TERM2GENE = pathwayToMetabolite, 
            TERM2NAME = pathwayToMetaboliteName)

head(summary(x))

metPathwayEnrichment<-as.data.frame(x)

}


#####
# run GSEA on 85 KEGG metabolic pathways
#####
library(dplyr)

filePath<-"~/work/Ed_lab/HCC_project/data/KEGG_metabolic_pathway/parsed_data"
fileName<-paste("KEGG_metabolic_pathway_list.Rd",sep="")
fileName<-file.path(filePath,fileName)

load(file=fileName)

if(FALSE){
fgseaRes_kegg_metabolic_pathway <- fgsea(pathways = metabolitePathway, 
                                         stats    = abs_signalToNoiseRatio_rank_KEGG_class_4,
                                         #stats    = signalToNoiseRatio_rank_class_4,
                                         eps      = 0.0,
                                         scoreType = "pos",
                                         minSize  = 5,
                                         maxSize  = 500)
}

fgseaRes_kegg_metabolic_pathway <- fgsea(pathways = metabolitePathway, 
                                         #stats    = abs_signalToNoiseRatio_rank_class_4,
                                         stats    = signalToNoiseRatio_rank_KEGG_class_1,
                                         eps      = 0.0,
                                         #scoreType = "pos",
                                         minSize  = 5,
                                         maxSize  = 500)

fgseaRes_kegg_metabolic_pathway<-fgseaRes_kegg_metabolic_pathway[order(padj, -abs(NES)), ]

head(fgseaRes_kegg_metabolic_pathway[order(padj, -abs(NES)), ], n=10)

########

# plotEnrichment(metabolitePathway[["Citrate cycle (TCA cycle)"]], signalToNoiseRatio_rank_KEGG_class_4)

######

fgseaResTidy <- fgseaRes_kegg_metabolic_pathway %>%
  as_tibble() %>%
  arrange( desc(abs(NES)) )

fgseaResTidy$rank<-seq(1,nrow(fgseaResTidy),1)

if(FALSE){
  # Show in a nice table:
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES) %>% 
    arrange(padj) %>% 
    DT::datatable()
}
#####

xAxisRange<-round(nrow(fgseaResTidy)/10)*10
xBreaks<-seq(10,xAxisRange,10)

graph <- ggplot(fgseaResTidy, aes(rank, NES, label=pathway))
graph <- graph + geom_point(size=1,col="grey50")
graph <- graph + geom_point(data=fgseaResTidy[fgseaResTidy$NES>0 & fgseaResTidy$padj<0.1, ],col="red", size=1)
graph <- graph + geom_point(data=fgseaResTidy[fgseaResTidy$NES<0 & fgseaResTidy$padj<0.1, ],col="royalblue", size=1)
graph <- graph + geom_hline(aes(yintercept=0))
graph <- graph + expand_limits(y=c(-3,3))
graph <- graph + scale_x_continuous(breaks=xBreaks)
graph <- graph + scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3))
#graph <- graph + scale_x_continuous(expand=c(0,0))
#graph <- graph + geom_col(aes(fill=padj<0.05))
#graph <- graph + coord_flip()
graph <- graph + geom_text_repel(data=fgseaResTidy[fgseaResTidy$padj<0.1,],
                                 #face="plain", family = 'ArialMT', color="black",
                                 size= 6 / .pt,
                                 min.segment.length = 0,
                                 box.padding = 0.5,
                                 max.overlaps = Inf,
                                 #nudge_x = 0.15,
                                 #nudge_y = 0.15,
                                 #segment.curvature = -0.1,
                                 #segment.ncp = 3,
                                 #segment.angle = 20
                                 force= 0.1,
                                 nudge_x      = 10,
                                 nudge_y = 1,
                                 #direction    = "y",
                                 direction = "y",
                                 hjust        = 0,
                                 #vjust = 0,
                                 segment.size = 0.3,
                                 segment.color = "black"
)
graph <- graph + labs(x="KEGG Metabolic Pathways", 
                      y="Normalized Enrichment Score",
                      title="class 1 vs others (metabolite)")
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=6, color="black"),
                       axis.text = element_text(size = 6, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=6, face="plain"),
                       axis.title.y = element_text(size=6, face="plain"),
                       axis.text.x = element_text(size=6, angle=0, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       #panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)


filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_metabolite_cc_cluster_GSEA_85_KEGG_metabolics_pathways_class_1_vs_others.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.25,width=2.25)
print(graph)
dev.off()

######
# GSEA for class 2 vs others
######

if(FALSE){
  fgseaRes_kegg_metabolic_pathway <- fgsea(pathways = metabolitePathway, 
                                           stats    = abs_signalToNoiseRatio_rank_KEGG_class_4,
                                           #stats    = signalToNoiseRatio_rank_class_4,
                                           eps      = 0.0,
                                           scoreType = "pos",
                                           minSize  = 5,
                                           maxSize  = 500)
}

fgseaRes_kegg_metabolic_pathway <- fgsea(pathways = metabolitePathway, 
                                         #stats    = abs_signalToNoiseRatio_rank_class_4,
                                         stats    = signalToNoiseRatio_rank_KEGG_class_2,
                                         eps      = 0.0,
                                         #scoreType = "pos",
                                         minSize  = 5,
                                         maxSize  = 500)

fgseaRes_kegg_metabolic_pathway<-fgseaRes_kegg_metabolic_pathway[order(padj, -abs(NES)), ]

head(fgseaRes_kegg_metabolic_pathway[order(padj, -abs(NES)), ], n=10)

########

# plotEnrichment(metabolitePathway[["Citrate cycle (TCA cycle)"]], signalToNoiseRatio_rank_KEGG_class_4)

######

fgseaResTidy <- fgseaRes_kegg_metabolic_pathway %>%
  as_tibble() %>%
  arrange( desc(abs(NES)) )

fgseaResTidy$rank<-seq(1,nrow(fgseaResTidy),1)

if(FALSE){
  # Show in a nice table:
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES) %>% 
    arrange(padj) %>% 
    DT::datatable()
}
#####

xAxisRange<-round(nrow(fgseaResTidy)/10)*10
xBreaks<-seq(10,xAxisRange,10)

graph <- ggplot(fgseaResTidy, aes(rank, NES, label=pathway))
graph <- graph + geom_point(size=1,col="grey50")
graph <- graph + geom_point(data=fgseaResTidy[fgseaResTidy$NES>0 & fgseaResTidy$padj<0.1, ],col="red", size=1)
graph <- graph + geom_point(data=fgseaResTidy[fgseaResTidy$NES<0 & fgseaResTidy$padj<0.1, ],col="royalblue", size=1)
graph <- graph + geom_hline(aes(yintercept=0))
graph <- graph + expand_limits(y=c(-3,3))
graph <- graph + scale_x_continuous(breaks=xBreaks)
graph <- graph + scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3))
#graph <- graph + scale_x_continuous(expand=c(0,0))
#graph <- graph + geom_col(aes(fill=padj<0.05))
#graph <- graph + coord_flip()
graph <- graph + geom_text_repel(data=fgseaResTidy[fgseaResTidy$padj<0.1,],
                                 #face="plain", family = 'ArialMT', color="black",
                                 size= 6 / .pt,
                                 min.segment.length = 0,
                                 box.padding = 0.5,
                                 max.overlaps = Inf,
                                 #nudge_x = 0.15,
                                 #nudge_y = 0.15,
                                 #segment.curvature = -0.1,
                                 #segment.ncp = 3,
                                 #segment.angle = 20
                                 force= 0.1,
                                 nudge_x      = 10,
                                 nudge_y = 1,
                                 #direction    = "y",
                                 direction = "y",
                                 hjust        = 0,
                                 #vjust = 0,
                                 segment.size = 0.3,
                                 segment.color = "black"
)
graph <- graph + labs(x="KEGG Metabolic Pathways", 
                      y="Normalized Enrichment Score",
                      title="class 2 vs others (metabolite)")
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=6, color="black"),
                       axis.text = element_text(size = 6, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=6, face="plain"),
                       axis.title.y = element_text(size=6, face="plain"),
                       axis.text.x = element_text(size=6, angle=0, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       #panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_metabolite_cc_cluster_GSEA_85_KEGG_metabolics_pathways_class_2_vs_others.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.25,width=2.25)
print(graph)
dev.off()


######
# run GSEA class 3 vs others
######

if(FALSE){
  fgseaRes_kegg_metabolic_pathway <- fgsea(pathways = metabolitePathway, 
                                           stats    = abs_signalToNoiseRatio_rank_KEGG_class_4,
                                           #stats    = signalToNoiseRatio_rank_class_4,
                                           eps      = 0.0,
                                           scoreType = "pos",
                                           minSize  = 5,
                                           maxSize  = 500)
}

fgseaRes_kegg_metabolic_pathway <- fgsea(pathways = metabolitePathway, 
                                         #stats    = abs_signalToNoiseRatio_rank_class_4,
                                         stats    = signalToNoiseRatio_rank_KEGG_class_3,
                                         eps      = 0.0,
                                         #scoreType = "pos",
                                         minSize  = 5,
                                         maxSize  = 500)

fgseaRes_kegg_metabolic_pathway<-fgseaRes_kegg_metabolic_pathway[order(padj, -abs(NES)), ]

head(fgseaRes_kegg_metabolic_pathway[order(padj, -abs(NES)), ], n=10)

########

# plotEnrichment(metabolitePathway[["Citrate cycle (TCA cycle)"]], signalToNoiseRatio_rank_KEGG_class_4)

######

fgseaResTidy <- fgseaRes_kegg_metabolic_pathway %>%
  as_tibble() %>%
  arrange( desc(abs(NES)) )

fgseaResTidy$rank<-seq(1,nrow(fgseaResTidy),1)

if(FALSE){
  # Show in a nice table:
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES) %>% 
    arrange(padj) %>% 
    DT::datatable()
}
#####

xAxisRange<-round(nrow(fgseaResTidy)/10)*10
xBreaks<-seq(10,xAxisRange,10)

graph <- ggplot(fgseaResTidy, aes(rank, NES, label=pathway))
graph <- graph + geom_point(size=1,col="grey50")
graph <- graph + geom_point(data=fgseaResTidy[fgseaResTidy$NES>0 & fgseaResTidy$padj<0.1, ],col="red", size=1)
graph <- graph + geom_point(data=fgseaResTidy[fgseaResTidy$NES<0 & fgseaResTidy$padj<0.1, ],col="royalblue", size=1)
graph <- graph + geom_hline(aes(yintercept=0))
graph <- graph + expand_limits(y=c(-3,3))
graph <- graph + scale_x_continuous(breaks=xBreaks)
graph <- graph + scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3))
#graph <- graph + scale_x_continuous(expand=c(0,0))
#graph <- graph + geom_col(aes(fill=padj<0.05))
#graph <- graph + coord_flip()
graph <- graph + geom_text_repel(data=fgseaResTidy[fgseaResTidy$padj<0.1,],
                                 #face="plain", family = 'ArialMT', color="black",
                                 size= 6 / .pt,
                                 min.segment.length = 0,
                                 box.padding = 0.5,
                                 max.overlaps = Inf,
                                 #nudge_x = 0.15,
                                 #nudge_y = 0.15,
                                 #segment.curvature = -0.1,
                                 #segment.ncp = 3,
                                 #segment.angle = 20
                                 force= 0.1,
                                 nudge_x      = 10,
                                 nudge_y = 1,
                                 #direction    = "y",
                                 direction = "y",
                                 hjust        = 0,
                                 #vjust = 0,
                                 segment.size = 0.3,
                                 segment.color = "black"
)
graph <- graph + labs(x="KEGG Metabolic Pathways", 
                      y="Normalized Enrichment Score",
                      title="class 3 vs others (metabolite)")
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=6, color="black"),
                       axis.text = element_text(size = 6, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=6, face="plain"),
                       axis.title.y = element_text(size=6, face="plain"),
                       axis.text.x = element_text(size=6, angle=0, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       #panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_metabolite_cc_cluster_GSEA_85_KEGG_metabolics_pathways_class_3_vs_others.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.25,width=2.25)
print(graph)
dev.off()


#####
# run GSEA class 4 vs others
#####

if(FALSE){
  fgseaRes_kegg_metabolic_pathway <- fgsea(pathways = metabolitePathway, 
                                           stats    = abs_signalToNoiseRatio_rank_KEGG_class_4,
                                           #stats    = signalToNoiseRatio_rank_class_4,
                                           eps      = 0.0,
                                           scoreType = "pos",
                                           minSize  = 5,
                                           maxSize  = 500)
}

fgseaRes_kegg_metabolic_pathway <- fgsea(pathways = metabolitePathway, 
                                         #stats    = abs_signalToNoiseRatio_rank_class_4,
                                         stats    = signalToNoiseRatio_rank_KEGG_class_4,
                                         eps      = 0.0,
                                         #scoreType = "pos",
                                         minSize  = 3,
                                         maxSize  = 500)

fgseaRes_kegg_metabolic_pathway<-fgseaRes_kegg_metabolic_pathway[order(padj, -abs(NES)), ]

head(fgseaRes_kegg_metabolic_pathway[order(padj, -abs(NES)), ], n=10)

########

# plotEnrichment(metabolitePathway[["Citrate cycle (TCA cycle)"]], signalToNoiseRatio_rank_KEGG_class_4)

######

fgseaResTidy <- fgseaRes_kegg_metabolic_pathway %>%
  as_tibble() %>%
  arrange( desc(abs(NES)) )

fgseaResTidy$rank<-seq(1,nrow(fgseaResTidy),1)

if(FALSE){
  # Show in a nice table:
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES) %>% 
    arrange(padj) %>% 
    DT::datatable()
}
#####

xAxisRange<-round(nrow(fgseaResTidy)/10)*10
xBreaks<-seq(10,xAxisRange,10)

graph <- ggplot(fgseaResTidy, aes(rank, NES, label=pathway))
graph <- graph + geom_point(size=1,col="grey50")
graph <- graph + geom_point(data=fgseaResTidy[fgseaResTidy$NES>0 & fgseaResTidy$padj<0.1, ],col="red", size=1)
graph <- graph + geom_point(data=fgseaResTidy[fgseaResTidy$NES<0 & fgseaResTidy$padj<0.1, ],col="royalblue", size=1)
graph <- graph + geom_hline(aes(yintercept=0))
graph <- graph + expand_limits(y=c(-3,3), x=c(0,xAxisRange))
graph <- graph + scale_x_continuous(breaks=xBreaks)
graph <- graph + scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3))
#graph <- graph + scale_x_continuous(expand=c(0,0))
#graph <- graph + geom_col(aes(fill=padj<0.05))
#graph <- graph + coord_flip()
graph <- graph + geom_text_repel(data=fgseaResTidy[fgseaResTidy$padj<0.1,],
                                 #face="plain", family = 'ArialMT', color="black",
                                 size= 6 / .pt,
                                 min.segment.length = 0,
                                 box.padding = 0.5,
                                 max.overlaps = Inf,
                                 #nudge_x = 0.15,
                                 #nudge_y = 0.15,
                                 #segment.curvature = -0.1,
                                 #segment.ncp = 3,
                                 #segment.angle = 20
                                 force= 0.1,
                                 nudge_x      = 10,
                                 nudge_y = 1,
                                 #direction    = "y",
                                 direction = "y",
                                 hjust        = 0,
                                 #vjust = 0,
                                 segment.size = 0.3,
                                 segment.color = "black"
)
graph <- graph + labs(x="KEGG Metabolic Pathways", 
                      y="Normalized Enrichment Score",
                      title="class 4 vs others (metabolite)")
graph <- graph + theme_classic()
graph <- graph + theme(text = element_text(size=6, color="black"),
                       axis.text = element_text(size = 6, color="black"),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=6, face="plain"),
                       axis.title.y = element_text(size=6, face="plain"),
                       axis.text.x = element_text(size=6, angle=0, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       #panel.grid.major = element_blank(), 
                       #panel.grid.minor = element_blank(),
                       #panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_line(colour = "grey",size=0.1),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/integrated_clustering/07_09_2021/combine_mRNA_met_cluster/metabolomics"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_metabolite_cc_cluster_GSEA_85_KEGG_metabolics_pathways_class_4_vs_others.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.25,width=2.25)
print(graph)
dev.off()





######
library(enrichplot)
library(clusterProfiler)

filePath<-"~/work/Ed_lab/HCC_project/data/KEGG_metabolic_pathway/parsed_data"
fileName<-paste("KEGG_metabolic_pathway_metaboliteInfoAll.txt",sep="")
#fileName<-paste("KEGG_metabolic_pathway_geneInfoAll.txt",sep="")
fileName<-file.path(filePath,fileName)

pathwayKEGG_metabolite<-fread(fileName,data.table=FALSE)
pathwayKEGG_metabolite<-pathwayKEGG_metabolite[,c("pathwayName","keggMetaboliteID")]
colnames(pathwayKEGG_metabolite)<-c("term","gene")

gseaRes_kegg_metabolic_pathway<-GSEA(signalToNoiseRatio_rank_KEGG_class_4,
                                     TERM2GENE = pathwayKEGG_metabolite,
                                     eps      = 0.0,
                                     pvalueCutoff = 1,
                                     pAdjustMethod = "BH",
                                     minGSSize = 5,
                                     maxGSSize = 500
)

gseaRes_kegg<-as.data.frame(gseaRes_kegg_metabolic_pathway)

dotplot(gseaRes_kegg_metabolic_pathway)

gseaplot2(gseaRes_kegg_metabolic_pathway,
          geneSetID=1, 
          title = gseaRes_kegg_metabolic_pathway$Description[1], 
          base_size=15)

gseaplot2(gseaRes_kegg_metabolic_pathway,
          geneSetID=2, 
          title = gseaRes_kegg_metabolic_pathway$Description[2], 
          base_size=15)
