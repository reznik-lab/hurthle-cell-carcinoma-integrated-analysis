#install.packages("devtools")
#devtools::install_github("cansysbio/ConsensusTME")
library(org.Hs.eg.db)
library(annotate)

######
filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/results/integrated_clustering/02_21_2021/combine_mRNA_met_cluster"

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
# load name mapping from histogy sample name to metabolomics sample name
filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/hcc_cancercell"
fileName<-"clin_to_YD_sampleName.csv"
fileName<-file.path(filePath,fileName)

nameMapping<-read.table(file=fileName,header=TRUE,sep=",",stringsAsFactors = FALSE)
nameMapping<-nameMapping[(nchar(nameMapping$name)!=0),]
colnames(nameMapping)<-c("Sample","Sample_HCC")

nameMapping<-nameMapping[grepl("T$",nameMapping$Sample),]


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


######
filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_sample_info_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

# load sampinfo
load(file=fileName)

sampleInfo<-data.frame(t(sampinfo),stringsAsFactors = FALSE)

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

#####

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/hcc_cancercell"
fileName<-"clin_to_YD_sampleName.csv"
fileName<-file.path(filePath,fileName)

nameMapping<-read.table(file=fileName,header=TRUE,sep=",",stringsAsFactors = FALSE)
colnames(nameMapping)<-c("Sample","Sample_HCC")

nameMapping_normal<-nameMapping[grepl("NL$",nameMapping$Sample),]
nameMapping<-nameMapping[(grepl("T$",nameMapping$Sample) & nchar(nameMapping$Sample_HCC)!=0),]

######
filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/hcc_cancercell"
fileName<-"HCC.tpm.csv"
#fileName<-"HCC.raw.csv"
fileName<-file.path(filePath,fileName)

data<-read.table(file=fileName,header=TRUE,sep=",",stringsAsFactors = FALSE)



#transcript_ID<-rownames(data)
transcript_ID<-data$Entrez_ID
rna_subset<-data[,nameMapping$Sample]
rownames(rna_subset)<-transcript_ID
colnames(rna_subset)<-nameMapping$Sample_HCC

rna_subset_normal<-data[,nameMapping_normal$Sample]
rownames(rna_subset_normal)<-transcript_ID
colnames(rna_subset_normal)<-nameMapping_normal$Sample


# include normal samples
data_subset<-cbind(rna_subset,rna_subset_normal)
#data_subset<-rna_subset


sdValue<-apply(data_subset,1,sd, na.rm=TRUE)
removedNames<-names(sdValue[ sdValue==0 | is.na(sdValue) ])
message(sprintf("Remove %s metbolite with no variation",length(removedNames)))
data_subset<-data_subset[!(rownames(data_subset) %in% removedNames),]
rna_subset<-data_subset
#rna_subset<-rna_subset[,ixsamples]

#######

# Change to hugo names
newnames = unlist( lookUp(rownames(rna_subset),'org.Hs.eg','SYMBOL') )
newnames2 = newnames[-which(is.na(newnames))]
newnames2_dupnames = which(duplicated(newnames2))

if (length(newnames2_dupnames) > 0)
{
  newnames2 = newnames2[-newnames2_dupnames,]
} # apparently unnecessary, no dups

dtrim = rna_subset[which(rownames(rna_subset) %in% names(newnames2)),]
rownames(dtrim) = newnames2[rownames(dtrim)]

rna_subset<-dtrim
#######


log2_rna_subset<-log2(rna_subset+1)

log2_rna_subset<-as.matrix(log2_rna_subset)
#######

sampleInfo_subset<-sampleInfo[sampleInfo$CLIENT_IDENTIFIER %in% nameMapping$Sample_HCC,]

#######

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/hcc_cancercell"
#fileName<-"MergedDeconvolution.HCC.csv"
fileName<-"MergedDeconvolution_add_charoentong_set_02_02_2021.HCC.csv"
fileName<-file.path(filePath,fileName)  

# Grab the RNA deconvolution data
decon = read.csv(fileName,header = TRUE,row.names = 1)

decon_tumor<-decon[nameMapping$Sample,]
rownames(decon_tumor)<-nameMapping$Sample_HCC

decon_normal<-decon[nameMapping_normal$Sample,]

decon<-rbind(decon_tumor,decon_normal)


#####

decon$log2CYT<-log2(decon$CYT)

log2CYT<-decon$log2CYT

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

####

decon_immune_subset<-decon[,c(1:97)]


#####

decon_immune_subset<-t(decon_immune_subset)



#####
d<-dist(t(decon_immune_subset),method="euclidean")
hc_result<-hclust(d,method="ward.D2")

plot(hc_result,cex=0.6,hang=-1)

d<-dist(t(decon_immune_subset),method="euclidean")
hc_result<-hclust(d,method="ward.D2")

plot(hc_result,cex=0.6,hang=-1)

sub_grp <- cutree(hc_result, k = 3)

plot(hc_result,cex=0.6,hang=-1)
rect.hclust(hc_result, k = 3, border = 2:5)

#####


library(factoextra)
library(NbClust)

df<-t(tme_signatures)

fviz_cluster(list(data = df, cluster = sub_grp))


# Elbow method
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

########

# https://rdrr.io/github/jhmadsen/ClustTools/man/SSE.html

#remotes::install_github("jhmadsen/ClustTools")
library("ClustTools")

df<-t(tme_signatures)

d<-dist(df,method="euclidean")
cluster.obj<-hclust(d,method="ward.D2")

plot(cluster.obj,cex=0.6,hang=-1)

#sub_grp <- cutree(hc_result, k = 4)

output<-list()

for(cluster_threshold in 1:10){
  
  ## Cut the hierarchical clustering tree
  chosen.clusters <- cutree(cluster.obj, k= cluster_threshold)
  
  ## Evaluate the clustering results with 'SSE' and 'SST'
  clusters.SSE <- SSE(df, chosen.clusters)$sumWithin
  clusters.SST <- SST(df)
  
  ## Calculate the r-squared for your cluster solution
  var <- 1-(clusters.SSE/clusters.SST)
  
  output[[cluster_threshold]]<-data.frame(cluster_threshold,var,stringsAsFactors = FALSE)
  
}

output<-rbind.fill(output)

plot(output$cluster_threshold,output$var,type="l")


#######


if(FALSE){
  
  library(pheatmap)
  
  data_subset<-tme_signatures[,consensus_sample_order]
  
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
  data_subset_norm <- t(apply(data_subset, 1, cal_z_score))
  
  data_subset_norm[data_subset_norm>2]<-2
  data_subset_norm[data_subset_norm<(-2)]<-(-2)
  
  #data_subset_norm<-tme_signatures[,consensus_sample_order]
  
  #####
  
  met_mRNA_cluster_label<-met_mRNA_cluster_label[sampleInfo_subset$CLIENT_IDENTIFIER,]
  
  my_histology_col<-data.frame(histology=sampleInfo_subset$HISTOLOGY2,stringsAsFactors = FALSE)
  rownames(my_histology_col)<-sampleInfo_subset$CLIENT_IDENTIFIER
  
  my_genetic_col<-data.frame(clingen_subset[,c(1:7)],stringsAsFactors = FALSE)
  my_genetic_col<-my_genetic_col[sampleInfo_subset$CLIENT_IDENTIFIER,]
  
  my_met_cluster_col<-data.frame(met_cluster=met_mRNA_cluster_label$met_cluster,stringsAsFactors = FALSE)
  rownames(my_met_cluster_col)<-rownames(met_mRNA_cluster_label)
  
  my_mRNA_cluster_col<-data.frame(mRNA_cluster=met_mRNA_cluster_label$mRNA_cluster,stringsAsFactors = FALSE)
  rownames(my_mRNA_cluster_col)<-rownames(met_mRNA_cluster_label)
  
  my_annotation<-cbind(my_histology_col,my_genetic_col,my_met_cluster_col,my_mRNA_cluster_col)
  
  #####
  
  my_consensus_cluster_col<-data.frame(consensus_cluster=met_mRNA_cluster_label$class,stringsAsFactors = FALSE)
  rownames(my_consensus_cluster_col)<-rownames(met_mRNA_cluster_label)
  
  my_annotation<-my_annotation[rownames(my_consensus_cluster_col),]
  
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
                             met_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1"),
                             mRNA_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1"),
                             #decon_cluster=c("1"="red","2"="royalblue","3"="yellow"),
                             #decon_cluster=c("1"="red","2"="royalblue","3"="yellow","4"="green","5"="pink"),
                             consensus_cluster=c("1"="red","2"="royalblue","3"="mediumpurple1","4"="pink1"))
  
  
  paletteLength <- 15
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(min(data_subset_norm), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(data_subset_norm)/paletteLength, max(data_subset_norm), length.out=floor(paletteLength/2)))
  
  
  if(FALSE){
    mygraph<-pheatmap(data_subset, 
                      #color =myColor, 
                      #breaks=myBreaks,  
                      annotation_col = my_sample_col,
                      cluster_rows=TRUE, 
                      cluster_cols=FALSE, 
                      scale="row", 
                      fontsize=7,
                      border_color = "white")
  }
  
  mygraph<-pheatmap(data_subset_norm, 
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
  
  filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/02_25_2021/decon/heatmap"
  dir.create(filePath,recursive = TRUE)
  fileName<-paste("HCC_decon_immuneSignatures_NES_scale_heatmap.pdf",sep="")
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_CID_2_vs_3_AQR_CR_with_sig_pvalue.pdf"
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_before.pdf"
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_late.pdf"
  
  fileName<-file.path(filePath,fileName)
  
  pdf(fileName,height=5,width=8)
  print(mygraph)
  dev.off()
  
  mygraph<-pheatmap(data_subset_norm, 
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
  
  filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/02_25_2021/decon/heatmap"
  dir.create(filePath,recursive = TRUE)
  fileName<-paste("HCC_decon_immuneSignatures_zscore_scale_heatmap.pdf",sep="")
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_CID_2_vs_3_AQR_CR_with_sig_pvalue.pdf"
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_before.pdf"
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_late.pdf"
  
  fileName<-file.path(filePath,fileName)
  
  pdf(fileName,height=5,width=8)
  print(mygraph)
  dev.off()
  
  
  pheatmap(data_subset_norm,
           cluster_rows=TRUE,
           cluster_cols = TRUE,
           scale="none")
  
  pheatmap(decon_tumor_NES,
           cluster_rows=TRUE,
           cluster_cols = TRUE,
           scale="none")
  
  
}


########
library(ggplot2)

#data_subset<-tme_signatures[,consensus_sample_order]
data_subset<-decon_immune_subset

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(data_subset, 1, cal_z_score))

mat<-data_subset_norm

data_subset<-mat
sdValue<-apply(data_subset,1,sd, na.rm=TRUE)
removedNames<-names(sdValue[ sdValue==0 | is.na(sdValue) ])
message(sprintf("Remove %s metbolite with no variation",length(removedNames)))
data_subset<-data_subset[!(rownames(data_subset) %in% removedNames),]
mat<-data_subset

tme_signatures_list<-rownames(mat)

test_result<-list()


for(idx in 1:length(tme_signatures_list)){
  
  picked_score<-tme_signatures_list[idx]
  
  message(sprintf("Work on %s",picked_score))
  
  vector_normal<-mat[picked_score,nameMapping_normal$Sample]
  
  vector_class_1<-mat[picked_score,sample_class_1]
  vector_class_2<-mat[picked_score,sample_class_2]
  vector_class_3<-mat[picked_score,sample_class_3]
  vector_class_4<-mat[picked_score,sample_class_4]
  
  vector_tumor<-c(vector_class_1,
                  vector_class_2,
                  vector_class_3,
                  vector_class_4)
  
  
  mean_normal<-mean(vector_normal)
  mean_class_1<-mean(vector_class_1)
  mean_class_2<-mean(vector_class_2)
  mean_class_3<-mean(vector_class_3)
  mean_class_4<-mean(vector_class_4)
  
  mean_tumor_no_class_1<-mean(c(vector_class_2,vector_class_3,vector_class_4))
  mean_tumor_no_class_2<-mean(c(vector_class_1,vector_class_3,vector_class_4))
  mean_tumor_no_class_3<-mean(c(vector_class_1,vector_class_2,vector_class_4))
  mean_tumor_no_class_4<-mean(c(vector_class_1,vector_class_2,vector_class_3))
  
  mean_tumor_class_3_4<-mean(c(vector_class_3,vector_class_4))
  mean_tumor_class_1_4<-mean(c(vector_class_1,vector_class_4))
  mean_tumor_class_1_3<-mean(c(vector_class_1,vector_class_3))
  
  
  mean_tumor<-mean(vector_tumor)
  
  #fc_class_1<-mean_class_1/mean_normal
  #fc_class_2<-mean_class_2/mean_normal
  #fc_class_3<-mean_class_3/mean_normal
  #fc_class_4<-mean_class_4/mean_normal
  
  wilcox_test_tumor_vs_normal<-wilcox.test(vector_tumor,vector_normal)
  
  wilcox_class_1_vs_normal<-wilcox.test(vector_class_1,vector_normal)
  wilcox_class_2_vs_normal<-wilcox.test(vector_class_2,vector_normal)
  wilcox_class_3_vs_normal<-wilcox.test(vector_class_3,vector_normal)
  wilcox_class_4_vs_normal<-wilcox.test(vector_class_4,vector_normal)
  
  wilcox_class_1_vs_others<-wilcox.test(vector_class_1,c(vector_class_2,vector_class_3,vector_class_4))
  wilcox_class_2_vs_others<-wilcox.test(vector_class_2,c(vector_class_1,vector_class_3,vector_class_4))
  wilcox_class_3_vs_others<-wilcox.test(vector_class_3,c(vector_class_1,vector_class_2,vector_class_4))
  wilcox_class_4_vs_others<-wilcox.test(vector_class_4,c(vector_class_1,vector_class_2,vector_class_3))
  
  
  wilcox_class_1_vs_class_3_4<-wilcox.test(vector_class_2,c(vector_class_3,vector_class_4))
  wilcox_class_3_vs_class_1_4<-wilcox.test(vector_class_3,c(vector_class_1,vector_class_4))
  wilcox_class_4_vs_class_1_3<-wilcox.test(vector_class_4,c(vector_class_1,vector_class_3))
  
  pvalue_tumor_vs_normal<-wilcox_test_tumor_vs_normal$p.value
  
  pvalue_class_1_vs_normal<-wilcox_class_1_vs_normal$p.value
  pvalue_class_2_vs_normal<-wilcox_class_2_vs_normal$p.value
  pvalue_class_3_vs_normal<-wilcox_class_3_vs_normal$p.value
  pvalue_class_4_vs_normal<-wilcox_class_4_vs_normal$p.value
  
  pvalue_class_1_vs_others<-wilcox_class_1_vs_others$p.value
  pvalue_class_2_vs_others<-wilcox_class_2_vs_others$p.value
  pvalue_class_3_vs_others<-wilcox_class_3_vs_others$p.value
  pvalue_class_4_vs_others<-wilcox_class_4_vs_others$p.value
  
  pvalue_class_1_vs_class_3_4<-wilcox_class_1_vs_class_3_4$p.value
  pvalue_class_3_vs_class_1_4<-wilcox_class_3_vs_class_1_4$p.value
  pvalue_class_4_vs_class_1_3<-wilcox_class_4_vs_class_1_3$p.value
  

  test_result[[idx]]<-data.frame(picked_score,
                                 mean_normal,
                                 mean_tumor,
                                 mean_class_1,
                                 mean_class_2,
                                 mean_class_3,
                                 mean_class_4,
                                 mean_tumor_no_class_1,
                                 mean_tumor_no_class_2,
                                 mean_tumor_no_class_3,
                                 mean_tumor_no_class_4,
                                 mean_tumor_class_3_4,
                                 mean_tumor_class_1_4,
                                 mean_tumor_class_1_3,
                                 pvalue_tumor_vs_normal,
                                 pvalue_class_1_vs_normal,
                                 pvalue_class_2_vs_normal,
                                 pvalue_class_3_vs_normal,
                                 pvalue_class_4_vs_normal,
                                 pvalue_class_1_vs_others,
                                 pvalue_class_2_vs_others,
                                 pvalue_class_3_vs_others,
                                 pvalue_class_4_vs_others,
                                 pvalue_class_1_vs_class_3_4,
                                 pvalue_class_3_vs_class_1_4,
                                 pvalue_class_4_vs_class_1_3,
                                 stringsAsFactors = FALSE)
  
  data_normal<-data.frame(name=picked_score,value=vector_normal,group="normal",stringsAsFactors = FALSE)
  data_tumor<-data.frame(name=picked_score,value=vector_tumor,group="tumor",stringsAsFactors = FALSE)
  data_class_1<-data.frame(name=picked_score,value=vector_class_1,group="tumor_class_1",stringsAsFactors = FALSE)
  data_class_2<-data.frame(name=picked_score,value=vector_class_2,group="tumor_class_2",stringsAsFactors = FALSE)
  data_class_3<-data.frame(name=picked_score,value=vector_class_3,group="tumor_class_3",stringsAsFactors = FALSE)
  data_class_4<-data.frame(name=picked_score,value=vector_class_4,group="tumor_class_4",stringsAsFactors = FALSE)
  
  plotData<-rbind(data_normal,
                  data_tumor,
                  data_class_1,
                  data_class_2,
                  data_class_3,
                  data_class_4)
  
  plotData$group<-factor(plotData$group,levels=c("normal","tumor","tumor_class_1","tumor_class_2","tumor_class_3","tumor_class_4"))
  
  yAxisBreaks<-seq(-3,3,1)
  
  graph <- ggplot(plotData,aes(x=group,y=value,fill=group))
  graph <- graph + geom_boxplot(outlier.shape = NA, lwd=0.2, fatten=1)
  graph <- graph + geom_jitter(position=position_jitter(0.2), size=0.2)
  graph <- graph + ggtitle(paste(picked_score,sep=""))
  graph <- graph + xlab("groups")
  graph <- graph + ylab("z-score")
  graph <- graph + scale_y_continuous(#expand = c(0,0),
                                  breaks=yAxisBreaks)
  graph <- graph + expand_limits(y=c(-3,3))
  #graph <- graph + scale_color_manual(values=c("grey","red","royalblue","mediumpurple1","pink1"))
  graph <- graph + scale_fill_manual(values=c("grey","indianred1","pink1","red","royalblue","mediumpurple1"))
  graph <- graph + theme_classic()
  graph <- graph + theme(axis.text = element_text(size = 7),
                 #axis.title = element_text(size = 7, face="bold"),
                 #axis.title.x = element_text(size=7),
                 axis.title.x = element_blank(),
                 axis.title.y = element_text(size=7),
                 text = element_text(size=7),
                 axis.text.x=element_text(colour="black",angle=45, hjust=1),
                 axis.text.y=element_text(colour="black"),
                 axis.ticks=element_line(colour="black"),
                 #axis.text.x = element_text(angle=45,hjust=1), 
                 #axis.title.x = element_blank(),
                 #axis.title.y = element_blank(),
                 #panel.border = element_rect(linetype = "solid", colour = "black"),
                 panel.border = element_blank(),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line=element_line(colour="black"),
                 plot.margin= margin(t=5, r=5, b=5, l=5, "pt"),
                 plot.title = element_text(lineheight=1.5,size=7,hjust=0.5),
                 plot.subtitle = element_text(lineheight=1.5,size=7,hjust=0.5),
                 legend.position="none")
  
  #print(graph)
  
  filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
  dir.create(filePath,recursive = TRUE)
  fileName<-paste("HCC_",picked_score,"_decon_immuneSignatures_zscore_scale_boxplot.pdf",sep="")
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_CID_2_vs_3_AQR_CR_with_sig_pvalue.pdf"
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_before.pdf"
  #fileName<-"AML_merged_normalized_heatmap_by_group_time_late.pdf"
  
  fileName<-file.path(filePath,fileName)
  
  pdf(fileName,height=2,width=2)
  print(graph)
  dev.off()
  
}

library(plyr)

test_result<-rbind.fill(test_result)

test_result_subset<-test_result[!grepl("^CB",test_result$picked_score),]


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


#test_result_subset<-test_result_subset[test_result_subset$picked_score %in% selected_TME_signaures,]

if(FALSE){
test_result_subset$padj_tumor_vs_normal<-p.adjust(test_result_subset$pvalue_tumor_vs_normal,method="BH")
test_result_subset$padj_class_1<-p.adjust(test_result_subset$pvalue_class_1,method="BH")
test_result_subset$padj_class_2<-p.adjust(test_result_subset$pvalue_class_2,method="BH")
test_result_subset$padj_class_3<-p.adjust(test_result_subset$pvalue_class_3,method="BH")
test_result_subset$padj_class_4<-p.adjust(test_result_subset$pvalue_class_4,method="BH")
}


filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-paste("HCC_decon_immuneSignatures_zscore_scale_wilcox_test.txt",sep="")
#fileName<-"AML_merged_normalized_heatmap_by_group_time_CID_2_vs_3_AQR_CR_with_sig_pvalue.pdf"
#fileName<-"AML_merged_normalized_heatmap_by_group_time_before.pdf"
#fileName<-"AML_merged_normalized_heatmap_by_group_time_late.pdf"

fileName<-file.path(filePath,fileName)

write.table(test_result_subset,file=fileName,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)


#####

source('~/work/Ed_lab/code_dev/metabolomics_pipeline/diff_abundance_test/makeVolcanoPlot_for_TME.R')

abs_log2 <- function(x){
  x[x==0] <- 1
  si <- sign(x)
  si * log2(si*x)
}

res<-data.frame(test_result_subset[,c("picked_score","mean_normal","mean_tumor","pvalue_tumor_vs_normal")])
res$log2FC<-res$mean_tumor - res$mean_normal
res$padj<-res$pvalue_tumor_vs_normal
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,
                         pvalueCol = 'padj',
                         threshold_pvalue = 0.05,
                         threshold_log2fc = 0,
                         show_legend_position = "none")

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_vs_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2,width=3)
print(plot)
dev.off()

######

res<-data.frame(test_result_subset[,c("picked_score","mean_normal","mean_tumor","pvalue_tumor_vs_normal")])

res$log2FC<-res$mean_tumor - res$mean_normal
res$padj<-p.adjust(res$pvalue_tumor_vs_normal,method="BH")
#res$padj<-res$pvalue_tumor_vs_normal
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 0)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_vs_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=6,width=6)
print(plot)
dev.off()



#####

res<-data.frame(test_result_subset[,c("picked_score","mean_normal","mean_class_1","pvalue_class_1_vs_normal")])
res$log2FC<-res$mean_class_1 - res$mean_normal
res$padj<-res$pvalue_class_1_vs_normal
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,
                         pvalueCol = 'padj',
                         threshold_pvalue = 0.05,
                         threshold_log2fc = 0,
                         title = "Volcano plot",
                         xlab = "mean z-score difference",
                         ylab = "-log10( raw p-value)",
                         savepath = NULL,
                         pointSize=0.7,
                         numOfCandidateLabels=30,
                         show_legend_position="none")

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_1_vs_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=2.2)
print(plot)
dev.off()

#####

removedTMEset<-c("B.Anti.apoptosis",
                 "B.Cytokines",
                 "B.Germinal.center",
                 "B.Naive_Memory",
                 "B.Pro.apoptosis",  
                 "B.Proliferation",
                 "GO_BCR_Signaling",
                 "Trm.BRCA",
                 "GO_TCR_Signaling")


res<-data.frame(test_result_subset[,c("picked_score","mean_normal","mean_class_2","pvalue_class_2_vs_normal")])
res$log2FC<-res$mean_class_2 - res$mean_normal
res$padj<-res$pvalue_class_2_vs_normal
#res$padj<-p.adjust(res$pvalue_class_2_vs_normal,method="BH")
rownames(res)<-res$picked_score

res_class_2_vs_normal<-res

plot<-makeVolcanoPlotTME(res=res,
                         pvalueCol = 'padj',
                         threshold_pvalue = 0.05,
                         threshold_log2fc = 0,
                         title = "Volcano plot",
                         xlab = "mean z-score difference",
                         ylab = "-log10( raw p-value)",
                         savepath = NULL,
                         pointSize=0.7,
                         numOfCandidateLabels=30,
                         show_legend_position="none")

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_2_vs_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=2.2)
print(plot)
dev.off()

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_2_vs_normal_volcano_plot_large_figure.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()


#####

res<-data.frame(test_result_subset[,c("picked_score","mean_normal","mean_class_3","pvalue_class_3_vs_normal")])
res$log2FC<-res$mean_class_3 - res$mean_normal
res$padj<-res$pvalue_class_3_vs_normal
#res$padj<-p.adjust(res$pvalue_class_3_vs_normal,method="BH")
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,
                         pvalueCol = 'padj',
                         threshold_pvalue = 0.05,
                         threshold_log2fc = 0,
                         title = "Volcano plot",
                         xlab = "mean z-score difference",
                         ylab = "-log10( raw p-value)",
                         savepath = NULL,
                         pointSize=0.7,
                         numOfCandidateLabels=30,
                         show_legend_position="none")

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_3_vs_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=2.2)
print(plot)
dev.off()

#####

res<-data.frame(test_result_subset[,c("picked_score","mean_normal","mean_class_4","pvalue_class_4_vs_normal")])
res$log2FC<-res$mean_class_4 - res$mean_normal
res$padj<-res$pvalue_class_4_vs_normal
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,
                         pvalueCol = 'padj',
                         threshold_pvalue = 0.05,
                         threshold_log2fc = 0,
                         title = "Volcano plot",
                         xlab = "mean z-score difference",
                         ylab = "-log10( raw p-value)",
                         savepath = NULL,
                         pointSize=0.7,
                         numOfCandidateLabels=30,
                         show_legend_position="none")

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_4_vs_normal_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=2.2)
print(plot)
dev.off()

#####

res<-data.frame(test_result_subset[,c("picked_score","mean_tumor_no_class_1","mean_class_1","pvalue_class_1_vs_others")])
res$log2FC<-res$mean_class_1 - res$mean_tumor_no_class_1
res$padj<-res$pvalue_class_1_vs_others
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 0)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_1_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=2.5,width=2.5)
print(plot)
dev.off()

#####

res<-data.frame(test_result_subset[,c("picked_score","mean_tumor_no_class_2","mean_class_2","pvalue_class_2_vs_others")])
res$log2FC<-res$mean_class_2 - res$mean_tumor_no_class_2
res$padj<-res$pvalue_class_2_vs_others
rownames(res)<-res$picked_score

res_class_2_vs_others<-res

plot<-makeVolcanoPlotTME(res=res,pvalueCol = 'padj',
                         threshold_pvalue = 0.05,
                         threshold_log2fc = 0, 
                         numOfCandidateLabels=5)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_2_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=3,width=3)
print(plot)
dev.off()

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_2_vs_others_volcano_plot_large_figure.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=10,width=10)
print(plot)
dev.off()


#####

res<-data.frame(test_result_subset[,c("picked_score","mean_tumor_no_class_3","mean_class_3","pvalue_class_3_vs_others")])
res$log2FC<-res$mean_class_3 - res$mean_tumor_no_class_3
res$padj<-res$pvalue_class_3_vs_others
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 0)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_3_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=6,width=6)
print(plot)
dev.off()

#####

res<-data.frame(test_result_subset[,c("picked_score","mean_tumor_no_class_4","mean_class_4","pvalue_class_4_vs_others")])
res$log2FC<-res$mean_class_4 - res$mean_tumor_no_class_4
res$padj<-res$pvalue_class_4_vs_others
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 0)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_4_vs_others_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=6,width=6)
print(plot)
dev.off()

#####

res<-data.frame(test_result_subset[,c("picked_score","mean_tumor_class_3_4","mean_class_1","pvalue_class_1_vs_class_3_4")])
res$log2FC<-res$mean_class_1 - res$mean_tumor_class_3_4
res$padj<-res$pvalue_class_1_vs_class_3_4
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 0)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_1_vs_class_3_4_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=6,width=6)
print(plot)
dev.off()

#####

res<-data.frame(test_result_subset[,c("picked_score","mean_tumor_class_1_4","mean_class_3","pvalue_class_3_vs_class_1_4")])
res$log2FC<-res$mean_class_3 - res$mean_tumor_class_1_4
res$padj<-res$pvalue_class_3_vs_class_1_4
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 0)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_3_vs_class_1_4_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=6,width=6)
print(plot)
dev.off()

#####

res<-data.frame(test_result_subset[,c("picked_score","mean_tumor_class_1_3","mean_class_4","pvalue_class_4_vs_class_1_3")])
res$log2FC<-res$mean_class_4 - res$mean_tumor_class_1_3
res$padj<-res$pvalue_class_4_vs_class_1_3
rownames(res)<-res$picked_score

plot<-makeVolcanoPlotTME(res=res,pvalueCol = 'padj',threshold_pvalue = 0.05,threshold_log2fc = 0)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_4_vs_class_1_3_volcano_plot.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=6,width=6)
print(plot)
dev.off()

#####
library(ggvenn)

sig_set<-list(class_2_vs_normal=res_class_2_vs_normal[res_class_2_vs_normal$pvalue_class_2_vs_normal<0.05,]$picked_score,
           class_2_vs_others=res_class_2_vs_others[res_class_2_vs_others$pvalue_class_2_vs_others<0.05,]$picked_score
           )

venn_graph<-ggvenn(sig_set, show_percentage = FALSE)

filePath<-"~/work/Ed_lab/HCC_project/results/TME_clustering/03_25_2021/decon/boxplot_zscore_scale"
dir.create(filePath,recursive = TRUE)
fileName<-"HCC_TME_integrated_clustering_tumor_class_2_venn_diagram.pdf"
fileName<-file.path(filePath,fileName)

pdf(fileName,height=6,width=6)
print(venn_graph)
dev.off()


#####
library(gplots)
# https://davetang.org/muse/2017/02/19/intersect-in-r/

my_venn<-venn(sig_set)
attr(x = my_venn, "intersections")



