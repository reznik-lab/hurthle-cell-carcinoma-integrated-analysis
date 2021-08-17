library(ggplot2)
#library(ggtern)
library(ggrepel)
library(data.table)
library(ggvenn)
library(plyr)
library(readxl)

#####
# load metabolite info
#####

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_metabolite_info_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

load(file=fileName)




#####

date<-"04_26_2021"

# p-value threshold
#threshold_pvalue = 0.1
threshold_pvalue = 0.05

# log2 fold change threshold
threshold_log2fc = 1

# Read in the various histologies
# DA: differential analysis

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_vs_normal.txt"
fileName<-file.path(filePath,fileName)
all_tumors_vs_normal = read.table(file=fileName,sep="\t",header = TRUE,quote="",fill=NA,stringsAsFactors = FALSE)
rownames(all_tumors_vs_normal)<-all_tumors_vs_normal$name


filePath<-"/Users/lium2/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/04_26_2021"
fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_tumor_vs_hcc_normal.txt"
fileName<-file.path(filePath,fileName)
hcc = read.table(file=fileName,sep="\t",header = TRUE,quote="",fill=NA,stringsAsFactors = FALSE)
rownames(hcc)<-hcc$name

fileName<-"HCC_differential_metabolite_abundance_wilcox_test_pd_ptctv_tumor_vs_pd_ptctv_normal.txt"
fileName<-file.path(filePath,fileName)
pd_and_ptctv = read.table(file=fileName,sep="\t",header = TRUE,quote="",fill=NA,stringsAsFactors = FALSE)
rownames(pd_and_ptctv)<-pd_and_ptctv$name

fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_vs_pd_and_ptctv.txt"
fileName<-file.path(filePath,fileName)
hcc_vs_pd_and_ptctv = read.table(file=fileName,sep="\t",header = TRUE,quote="",fill=NA,stringsAsFactors = FALSE)
rownames(hcc_vs_pd_and_ptctv)<-hcc_vs_pd_and_ptctv$name

#####
# load supplementary table 1 from Xu et al J Proteome Res 2015
# PMID: 26130307
#####

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/Xu_et_al_J_Proteome_Res_2015"
fileName<-"parsed_supplementary_tables.xlsx"
fileName<-file.path(filePath,fileName)
ptc_vs_normal = read_excel(path=fileName,sheet="supplementary table 1")
ptc_vs_normal<-as.data.frame(ptc_vs_normal)
rownames(ptc_vs_normal)<-ptc_vs_normal$`Metabolite Name`

ptc_vs_normal$log2FC<-log2(ptc_vs_normal$FC)
#ptc_vs_normal$name<-tolower(ptc_vs_normal$`Metabolite Name`)

# 35 metabolites overlaps
overlap_KEGG<-intersect(ptc_vs_normal$KEGG,metinfo$H_KEGG)

# 11 metabolites non-overlap
#setdiff(ptc_vs_normal$KEGG,metinfo$H_KEGG)

#ptc_vs_normal$H_name<-tolower(ptc_vs_normal$`Metabolite Name`)

metinfo_subset<-metinfo[,c("H_name","H_KEGG")]
colnames(metinfo_subset)<-c("H_name","KEGG")

ptc_vs_normal<-merge(ptc_vs_normal,metinfo_subset,by="KEGG",all.x = TRUE)

ptc_vs_normal$name<-ptc_vs_normal$H_name

#######
# handle differential metabolite abundance list
#######

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_metabolite_info_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

load(file=fileName)

metinfoSubset<-metinfo[,c("name","H_KEGG")]



#######
      
#merged_dat<-merge(hcc,hcc_vs_pd_and_ptctv,by="name",all.x = TRUE)
merged_dat<-merge(hcc,pd_and_ptctv,by="name",all.x = TRUE)
plot(merged_dat$log2FC.x,merged_dat$log2FC.y, pch=19)


hcc_common_pd_and_ptctv<-merged_dat[,c("name","log2FC.x","log2FC.y","pathway.x","padj.x","padj.y")]
colnames(hcc_common_pd_and_ptctv)<-c("name","log2FC_HCC","log2FC_other","pathway","padj_HCC","padj_other")      


merged_dat<-merge(hcc,hcc_vs_pd_and_ptctv,by="name",all.x = TRUE)
plot(merged_dat$log2FC.x,merged_dat$log2FC.y, pch=19)


hcc_specific_pd_and_ptctv<-merged_dat[,c("name","log2FC.x","log2FC.y","pathway.x","padj.x","padj.y")]
colnames(hcc_specific_pd_and_ptctv)<-c("name","log2FC_HCC","log2FC_other","pathway","padj_HCC","padj_other")      


######

hcc_common_pd_and_ptctv<-hcc_common_pd_and_ptctv[hcc_common_pd_and_ptctv$padj_HCC<0.05 & hcc_common_pd_and_ptctv$padj_other<0.05,]

hcc_specific_pd_and_ptctv<-hcc_specific_pd_and_ptctv[hcc_specific_pd_and_ptctv$padj_HCC<0.05 & hcc_specific_pd_and_ptctv$padj_other<0.05,]

hcc_to_normal<-hcc[hcc$pvalue<0.05,]

#metaboliteNameSet<-list("HCC_to_normal" = hcc_to_normal$name,
#                        "HCC_common" = hcc_common_pd_and_ptctv$name,
#                        "HCC_specific" = hcc_specific_pd_and_ptctv$name)


metaboliteNameSet<-list("HCC_vs_normal" = hcc_to_normal$name,
                        "pd_ptctv_vs_normal" = pd_and_ptctv[pd_and_ptctv$padj<0.05,]$name,
                        "HCC_vs_pd_ptctv" = hcc_vs_pd_and_ptctv[hcc_vs_pd_and_ptctv$padj<0.05,]$name)


graph<-ggvenn(metaboliteNameSet, show_percentage = FALSE)
print(graph)

filePath<-file.path("/Users/lium2/work/Ed_lab/HCC_project/results/comparative_metabolomics",date,"weighted_differential_abundance_score")
dir.create(filePath, showWarnings = TRUE, recursive = TRUE)
#fileName<-"HCC_differential_abundance_score.pdf"
fileName<-paste("HCC_sig_metabolites_three_set_venn_diagram.pdf",sep="")
#fileName<-"HCC_differential_abundance_score_show_less_than_5_measured_metabolites_in_pathway.pdf"

fileName<-file.path(filePath,fileName)

pdf(fileName,height=7,width=7)
print(graph)
dev.off()


library(gplots)

venn.table<-venn(metaboliteNameSet)

print(venn.table)

########


metaboliteNameSet<-list("HCC_vs_normal" = hcc_to_normal$name,
                        "pd_ptctv_vs_normal" = pd_and_ptctv[pd_and_ptctv$padj<0.05,]$name)

venn.table<-venn(metaboliteNameSet)

graph<-ggvenn(metaboliteNameSet, show_percentage = FALSE)

print(graph)

filePath<-file.path("/Users/lium2/work/Ed_lab/HCC_project/results/comparative_metabolomics",date,"weighted_differential_abundance_score")
dir.create(filePath, showWarnings = TRUE, recursive = TRUE)
#fileName<-"HCC_differential_abundance_score.pdf"
fileName<-paste("HCC_sig_metabolites_two_set_venn_diagram.pdf",sep="")
#fileName<-"HCC_differential_abundance_score_show_less_than_5_measured_metabolites_in_pathway.pdf"

fileName<-file.path(filePath,fileName)

pdf(fileName,height=7,width=7)
print(graph)
dev.off()


HCC_sig_metabolites_common<-attributes(venn.table)$intersections$`HCC_vs_normal:pd_ptctv_vs_normal`
HCC_sig_metabolites_specific<-attributes(venn.table)$intersections$HCC_vs_normal


matchedDiffMetAbundanceResult<-merge(hcc,metinfoSubset,by="name",all.x=TRUE)





########





########

plotData2<-hcc_common_pd_and_ptctv

#plotData2<-plotData2[abs(plotData2$log2FC_HCC)>1,]
#plotData2<-plotData2[abs(plotData2$log2FC_other)>1,]

#plotData3<-plotData2[!(plotData2$pathway %in% "Lipid"),]

criteria1<-abs(plotData2$log2FC_HCC)>4
criteria2<-plotData2$padj_HCC<0.05
criteria3<-abs(plotData2$log2FC_other)>4
criteria4<-plotData2$padj_other<0.05  
  
plotData3_1<-plotData2[criteria1 & criteria2 & criteria3 & criteria4 & plotData2$log2FC_HCC>4,]
plotData3_2<-plotData2[criteria1 & criteria2 & criteria3 & criteria4 & plotData2$log2FC_HCC< (-4),]


plotData3<-plotData2[criteria1 & criteria2 & criteria3 & criteria4, ]


graph <- ggplot(data=plotData2,aes(x=log2FC_HCC,y=log2FC_other, label=name))
graph <- graph + geom_point(color="grey", size = 1)
#graph <- graph + geom_point(data=plotData[plotData$padj_HCC<0.05 & plotData$padj_other<0.05,], color="red")
graph <- graph + geom_point(data=plotData3_1, 
                            size = 1,
                            color="red"
                            #color=pathway)
                            )

graph <- graph + geom_point(data=plotData3_2, 
                            size = 1,
                            color="royalblue"
                            #color=pathway)
)

graph <- graph + geom_vline(xintercept=0, linetype = "dashed", color="grey65")
graph <- graph + geom_hline(yintercept=0, linetype = "dashed", color ="grey65")
#graph <- graph + geom_abline(intercept = 0, slope =1)
graph <- graph + xlab("log2 fold change (HCC tumor/normal)")
graph <- graph + ylab("log2 fold change (PD + PTCTV/normal)")
graph <- graph + expand_limits(x=c(-11,11),y=c(-11,11))
graph <- graph + geom_text_repel(
                                  data=plotData3,
                                  size = 7 / .pt,
                                  color="black",
                                  min.segment.length = 0,
                                  segment.color = "grey50",
                                  #box.padding = 0.2,
                                  max.overlaps = Inf,
                                  segment.size = 0.2,
                                  #nudge_x = 0.15,
                                  #nudge_y = 0.15,
                                  #segment.curvature = -0.1,
                                  #segment.ncp = 3,
                                  #segment.angle = 20
                                )  
graph <- graph + theme_classic()
graph <- graph + theme(axis.text = element_text(size = 7),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=7),
                       axis.title.y = element_text(size=7),
                       axis.text.x = element_text(size=7, angle=0, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       #panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)

filePath<-"~/work/Ed_lab/HCC_project/results/comparative_metabolomics/03_10_2021/relative_to_normal"
dir.create(filePath,recursive = TRUE)
fileName<-paste("HCC_tumor_vs_pd_ptctv_to_normal_sig_metabolites.pdf",sep="")
fileName<-file.path(filePath,fileName)

pdf(fileName,width=3,height=3)
print(graph)
dev.off()

########

########

plotData2<-hcc_specific_pd_and_ptctv

#plotData2<-plotData2[abs(plotData2$log2FC_HCC)>1,]
#plotData2<-plotData2[abs(plotData2$log2FC_other)>1,]

#plotData3<-plotData2[!(plotData2$pathway %in% "Lipid"),]

criteria1<-abs(plotData2$log2FC_HCC)>2
criteria2<-plotData2$padj_HCC<0.05
criteria3<-abs(plotData2$log2FC_other)>2
criteria4<-plotData2$padj_other<0.05  

plotData3_1<-plotData2[criteria1 & criteria2 & criteria3 & criteria4 & plotData2$log2FC_HCC>2,]
plotData3_2<-plotData2[criteria1 & criteria2 & criteria3 & criteria4 & plotData2$log2FC_HCC< (-2),]


plotData3<-plotData2[criteria1 & criteria2 & criteria3 & criteria4, ]


graph <- ggplot(data=plotData2,aes(x=log2FC_HCC,y=log2FC_other, label=name))
graph <- graph + geom_point(color="grey", size = 1)
#graph <- graph + geom_point(data=plotData[plotData$padj_HCC<0.05 & plotData$padj_other<0.05,], color="red")
graph <- graph + geom_point(data=plotData3_1, 
                            size = 1,
                            color="red"
                            #color=pathway)
)

graph <- graph + geom_point(data=plotData3_2, 
                            size = 1,
                            color="royalblue"
                            #color=pathway)
)

graph <- graph + geom_vline(xintercept=0, linetype = "dashed", color="grey65")
graph <- graph + geom_hline(yintercept=0, linetype = "dashed", color ="grey65")
#graph <- graph + geom_abline(intercept = 0, slope =1)
graph <- graph + xlab("log2 fold change (HCC tumor/normal)")
graph <- graph + ylab("log2 fold change (HCC tumor/PD + PTCTV)")
graph <- graph + expand_limits(x=c(-11,11),y=c(-11,11))
graph <- graph + geom_text_repel(
  data=plotData3,
  size = 7 / .pt,
  color="black",
  min.segment.length = 0,
  segment.color = "grey50",
  #box.padding = 0.2,
  max.overlaps = Inf,
  segment.size = 0.2,
  #nudge_x = 0.15,
  #nudge_y = 0.15,
  #segment.curvature = -0.1,
  #segment.ncp = 3,
  #segment.angle = 20
)  
graph <- graph + theme_classic()
graph <- graph + theme(axis.text = element_text(size = 7),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=7),
                       axis.title.y = element_text(size=7),
                       axis.text.x = element_text(size=7, angle=0, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       #panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)

filePath<-"~/work/Ed_lab/HCC_project/results/comparative_metabolomics/04_26_2021/relative_to_normal"
dir.create(filePath,recursive = TRUE)
fileName<-paste("HCC_tumor_to_normal_vs_HCC_to_pd_ptctv_sig_metabolites.pdf",sep="")
fileName<-file.path(filePath,fileName)

pdf(fileName,width=3,height=3)
print(graph)
dev.off()

########
#  plot all tumor/normal to PTC/normal
########

library(dplyr)

merged_dat<-merge(all_tumors_vs_normal,ptc_vs_normal,by="name")
plot(merged_dat$log2FC.x,merged_dat$log2FC.y, pch=19)


all_tumors_and_ptc<-merged_dat[,c("name","log2FC.x","log2FC.y","pathway","pvalue","padj","P")]
colnames(all_tumors_and_ptc)<-c("name","log2FC_tumor","log2FC_ptc","pathway","pval_tumor","padj_tumor","pval_ptc")   




plotData2<-all_tumors_and_ptc


##########

all_tumor_up<-nrow(all_tumors_and_ptc[all_tumors_and_ptc$log2FC_tumor>0 & all_tumors_and_ptc$padj_tumor<0.05,])
all_tumor_down<-nrow(all_tumors_and_ptc[all_tumors_and_ptc$log2FC_tumor<0 & all_tumors_and_ptc$padj_tumor<0.05,])

ptc_up<-nrow(all_tumors_and_ptc[all_tumors_and_ptc$log2FC_ptc>0 & all_tumors_and_ptc$pval_ptc<0.05,])
ptc_down<-nrow(all_tumors_and_ptc[all_tumors_and_ptc$log2FC_ptc<0 & all_tumors_and_ptc$pval_ptc<0.05,])

contingencyTable<-data.frame(matrix(c(all_tumor_up,all_tumor_down,ptc_up,ptc_down),nrow=2))
colnames(contingencyTable)<-c("all_tumor","ptc")
rownames(contingencyTable)<-c("up","down")

fisher.test(contingencyTable,alternative = "two.sided")

#####
library(gplots)
df<-as.table(as.matrix(contingencyTable))

balloonplot(t(df))

chisq<-chisq.test(df)

#########


#plotData2<-plotData2[abs(plotData2$log2FC_HCC)>1,]
#plotData2<-plotData2[abs(plotData2$log2FC_other)>1,]

#plotData3<-plotData2[!(plotData2$pathway %in% "Lipid"),]

#criteria1<-abs(plotData2$log2FC_HCC)>2
#criteria2<-plotData2$padj_HCC<0.05
#criteria3<-abs(plotData2$log2FC_other)>2
#criteria4<-plotData2$padj_other<0.05  

#plotData3_1<-plotData2[criteria1 & criteria2 & criteria3 & criteria4 & plotData2$log2FC_HCC>2,]
#plotData3_2<-plotData2[criteria1 & criteria2 & criteria3 & criteria4 & plotData2$log2FC_HCC< (-2),]


#plotData3<-plotData2[criteria1 & criteria2 & criteria3 & criteria4, ]


criteria1<-plotData2$padj_tumor<0.05

plotData3<-plotData2[criteria1,]

graph <- ggplot(data=plotData2,aes(x=log2FC_tumor,y=log2FC_ptc, label=name))
graph <- graph + geom_point(color="grey", size = 1)
#graph <- graph + geom_point(data=plotData[plotData$padj_HCC<0.05 & plotData$padj_other<0.05,], color="red")
graph <- graph + geom_point(data=plotData3, 
                            size = 1,
                            color="black"
                            #color=pathway)
)

graph <- graph + geom_point(data=plotData3_2, 
                            size = 1,
                            color="royalblue"
                            #color=pathway)
)

graph <- graph + geom_vline(xintercept=0, linetype = "dashed", color="grey65")
graph <- graph + geom_hline(yintercept=0, linetype = "dashed", color ="grey65")
#graph <- graph + geom_abline(intercept = 0, slope =1)
graph <- graph + xlab("log2 fold change (HCC, HA, PD, PTCTV/normal)")
graph <- graph + ylab("log2 fold change (PTC/normal)")
graph <- graph + expand_limits(x=c(-5,5),y=c(-5,5))

graph <- graph + geom_text_repel(
  data=plotData2,
  size = 7 / .pt,
  color="black",
  min.segment.length = 0,
  segment.color = "grey50",
  #box.padding = 0.2,
  max.overlaps = Inf,
  segment.size = 0.2,
  #nudge_x = 0.15,
  #nudge_y = 0.15,
  #segment.curvature = -0.1,
  #segment.ncp = 3,
  #segment.angle = 20
)  
graph <- graph + theme_classic()
graph <- graph + theme(axis.text = element_text(size = 7),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=7),
                       axis.title.y = element_text(size=7),
                       axis.text.x = element_text(size=7, angle=0, hjust=1),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       #panel.grid.major = element_line(colour = "grey",size=0.1), 
                       #panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                       #legend.position="none")
                       legend.position = "right"
)

print(graph)

filePath<-"~/work/Ed_lab/HCC_project/results/comparative_metabolomics/06_22_2021/relative_to_normal"
dir.create(filePath,recursive = TRUE)
#fileName<-paste("all_tumor_to_normal_vs_PTC_to_normal_overlap_metabolites.pdf",sep="")
fileName<-paste("all_tumor_to_normal_vs_PTC_to_normal_overlap_metabolites_v2.pdf",sep="")
fileName<-file.path(filePath,fileName)

pdf(fileName,width=3,height=3)
print(graph)
dev.off()
