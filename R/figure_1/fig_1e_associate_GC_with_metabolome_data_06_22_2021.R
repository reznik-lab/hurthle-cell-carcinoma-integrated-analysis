library(ggplot2)
library(ggbeeswarm)
library(plyr)
library(ggrepel)
#library(xlsx)
library(readxl)

######


filePath<-"~/work/Ed_lab/HCC_project/data/From_Andy/07_28_2020"
fileName<-"mapping_with_metabolone_data.xlsx"  
fileName<-file.path(filePath,fileName)

sheet2use = 'Sheet1'
#tubeSampleInfo = read.xlsx(fileName,sheetName=sheet2use, check.names=FALSE,startRows=2, stringsAsFactors=FALSE)
#gc_metabolites = read.xlsx(fileName,sheetName=sheet2use, check.names=FALSE,stringsAsFactors=FALSE)
gc_metabolites = read_excel(fileName,sheet=sheet2use)
gc_metabolites<-as.data.frame(gc_metabolites)
rownames(gc_metabolites)<-gc_metabolites$name

metabolite_picked<-na.omit(gc_metabolites$H_name)

#######

filePath<-"/Users/lium2/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_met_Elisa_PNQ_normalized_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

# load met
load(file=fileName)

metabolon_met_log2<-log2(met)



#######

filePath<-"~/work/Ed_lab/HCC_project/data/From_Andy/07_28_2020"
fileName<-"GC data HCC samples_Jinsung_20200715.xlsx"  
fileName<-file.path(filePath,fileName)

sheet2use = 'combined data'
#tubeSampleInfo = read.xlsx(fileName,sheetName=sheet2use, check.names=FALSE,startRows=2, stringsAsFactors=FALSE)
#data_sheet1 = read.xlsx(fileName,sheetName=sheet2use, check.names=FALSE,stringsAsFactors=FALSE)
data_sheet1 = read_excel(fileName,sheet=sheet2use)
data_sheet1<-as.data.frame(data_sheet1)

sampleName<-c("YD4","YD5",
              "YD7","YD8",
              "YD12","YD13",
              "YD22","YD23",
              "YD24","YD25",
              "YD39","YD40",
              "YD54","YD55",
              "YD58","YD59",
              "YD60","YD61",
              "YD63","YD64",
              "YD66",
              "YD67","YD68",
              "YD70","YD71")

sampleName_tumor<-c("YD4",
              "YD7",
              "YD12",
              "YD22",
              "YD24",
              "YD39",
              "YD54",
              "YD58",
              "YD60",
              "YD63",
              "YD67",
              "YD70")

sampleName_normal<-c("YD5",
              "YD8",
              "YD13",
              "YD23",
              "YD25",
              "YD40",
              "YD55",
              "YD59",
              "YD61",
              "YD64",
              "YD66",
              "YD68",
              "YD71")




data_raws<-c(2:26)
sampleName_col<-c(3)

gc_metabolite<-list()

gc_metabolite[["pyruvate"]]<-data_sheet1[data_raws,c(sampleName_col,9)]

gc_metabolite[["lactate"]]<-data_sheet1[data_raws,c(sampleName_col,11)]

gc_metabolite[["alanine"]]<-data_sheet1[data_raws,c(sampleName_col,13)]
gc_metabolite[["valine"]]<-data_sheet1[data_raws,c(sampleName_col,15)]
gc_metabolite[["leucine"]]<-data_sheet1[data_raws,c(sampleName_col,17)]
gc_metabolite[["isoleucine"]]<-data_sheet1[data_raws,c(sampleName_col,19)]
gc_metabolite[["proline"]]<-data_sheet1[data_raws,c(sampleName_col,21)]
gc_metabolite[["glycine"]]<-data_sheet1[data_raws,c(sampleName_col,23)]
gc_metabolite[["succinate"]]<-data_sheet1[data_raws,c(sampleName_col,25)]
gc_metabolite[["fumarate"]]<-data_sheet1[data_raws,c(sampleName_col,27)]
gc_metabolite[["serine"]]<-data_sheet1[data_raws,c(sampleName_col,29)]
gc_metabolite[["threonine"]]<-data_sheet1[data_raws,c(sampleName_col,31)]
gc_metabolite[["malate"]]<-data_sheet1[data_raws,c(sampleName_col,33)]
gc_metabolite[["methionine"]]<-data_sheet1[data_raws,c(sampleName_col,35)]
gc_metabolite[["aspartate"]]<-data_sheet1[data_raws,c(sampleName_col,37)]
# what are cystein 322 and 220?
#cysteine 322
#cysteine 220
gc_metabolite[["d5_2HG"]]<-data_sheet1[data_raws,c(sampleName_col,43)]
gc_metabolite[["alpha-ketoglutarate"]]<-data_sheet1[data_raws,c(sampleName_col,45)]
gc_metabolite[["2-hydroxyglutarate"]]<-data_sheet1[data_raws,c(sampleName_col,47)]
gc_metabolite[["arginine"]]<-data_sheet1[data_raws,c(sampleName_col,49)]
gc_metabolite[["glutamate"]]<-data_sheet1[data_raws,c(sampleName_col,51)]
gc_metabolite[["phenylalanine"]]<-data_sheet1[data_raws,c(sampleName_col,53)]
gc_metabolite[["asparagine"]]<-data_sheet1[data_raws,c(sampleName_col,55)]
gc_metabolite[["glutamine"]]<-data_sheet1[data_raws,c(sampleName_col,57)]
gc_metabolite[["citrate"]]<-data_sheet1[data_raws,c(sampleName_col,59)]
gc_metabolite[["lysine"]]<-data_sheet1[data_raws,c(sampleName_col,61)]
gc_metabolite[["histidine"]]<-data_sheet1[data_raws,c(sampleName_col,63)]
gc_metabolite[["tyrosine"]]<-data_sheet1[data_raws,c(sampleName_col,65)]
gc_metabolite[["tryptophan"]]<-data_sheet1[data_raws,c(sampleName_col,67)]

# GC name / metabolon name
#gc_cis-aconitate / aconitate [cis or trans]
#ascorbate / ascorbate (Vitamin C)
#DL-2-aminoadipic acid / 2-aminoadipate
#3-hydroxyglutarate

#2-hydroxyadipate (no GC data)
#D-2-phosphoglyceric acid / 2-phosphoglycerate

# additional data
gc_metabolite[["aconitate [cis or trans]"]]<-data_sheet1[data_raws,c(sampleName_col,69)]
gc_metabolite[["ascorbate (Vitamin C)"]]<-data_sheet1[data_raws,c(sampleName_col,71)]
gc_metabolite[["2-aminoadipate"]]<-data_sheet1[data_raws,c(sampleName_col,73)]
gc_metabolite[["3-hydroxyglutarate"]]<-data_sheet1[data_raws,c(sampleName_col,75)]
gc_metabolite[["2-phosphoglycerate"]]<-data_sheet1[data_raws,c(sampleName_col,77)]

######

metabolon_met_log2_subset<-metabolon_met_log2[metabolite_picked,sampleName]

######

metaboliteList<-c("pyruvate",
                  "lactate",
                  "alanine",
                  "valine",
                  "leucine",
                  "isoleucine",
                  "proline",
                  "glycine",
                  "succinate",
                  "fumarate",
                  "serine",
                  "threonine",
                  "malate",
                  "methionine",
                  "aspartate",
                  #cysteine 322
                  #cysteine 220
                  #"d5_2HG",
                  "alpha-ketoglutarate",
                  "2-hydroxyglutarate",
                  "arginine",
                  "glutamate",
                  "phenylalanine",
                  "asparagine",
                  "glutamine",
                  "citrate",
                  "lysine",
                  "histidine",
                  "tyrosine",
                  "tryptophan",
                  "aconitate [cis or trans]",
                  "ascorbate (Vitamin C)",
                  "2-aminoadipate",
                  "3-hydroxyglutarate",
                  #"2-hydroxyadipate",
                  "2-phosphoglycerate"
                  )

res_summary<-list()

filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/tumor_and_normal"
dir.create(filePath,recursive = TRUE)

for(metabolite in metaboliteList){

#metabolite<-"pyruvate"

dat1<-gc_metabolite[[metabolite]]
dat1[dat1==0]<-NA
colnames(dat1)<-c("sampleName","gc_abundance")
dat1$sampleName<-do.call(rbind,strsplit(dat1$sampleName,"_"))[,1]
dat1$gc_abundance<-log2(as.numeric(dat1$gc_abundance))

met_vector<-metabolon_met_log2_subset[metabolite,]
dat2<-data.frame(names(met_vector),met_vector,stringsAsFactors = FALSE)
colnames(dat2)<-c("sampleName","metabolon_abundance")

cc<-merge(dat1,dat2,by="sampleName")
rownames(cc)<-cc$sampleName
cc<-cc[sampleName,]

maxValue<-round(max(cc$gc_abundance,cc$metabolon_abundance,na.rm=TRUE))+1
minValue<-round(min(cc$gc_abundance,cc$metabolon_abundance,na.rm=TRUE))-1


plotDat<-cc

#######
#library(ggpubr)

#ggscatter(plotDat,x="metabolon_abundance",y="gc_abundance",
#          add = "reg.line", conf.int = TRUE, 
#          cor.coef = TRUE, cor.method = "pearson",
#          xlab = "metabolon, log2(abundance)", ylab = "GC, log2(abundance)")
          
#######

res<-cor.test(plotDat$metabolon_abundance,plotDat$gc_abundance,method="pearson")
pvalue<-res$p.value
pvalue_for_plot<-formatC(res$p.value,format="e", digits=2)
coeff<-as.numeric(round(res$estimate,digit=2))

res_summary[[metabolite]]<-data.frame(metabolite,coeff,pvalue,stringsAsFactors = FALSE)


#######

#plot(log2(as.numeric(gc_pyruvate[,2])),metabolon_met_log2_subset["pyruvate",sampleName],pch=19,xlim=c(minValue,maxValue),ylim=c(minValue,maxValue))


graph <- ggplot(plotDat,aes(x=metabolon_abundance,y=gc_abundance))
graph <- graph + geom_point(size=1,shape=19)
graph <- graph + labs(title=metabolite)
graph <- graph + labs(subtitle=paste("R = ",coeff,", ","p = ",pvalue_for_plot,sep=""))
graph <- graph + geom_smooth(method=lm)
graph <- graph + xlab("log2(LC-MS abundance)")
graph <- graph + ylab("log2(GC-MS abundance)")

#graph <- graph + geom_abline(intercept = 0,slope=1,color="grey")
#graph <- graph + xlim(minValue,maxValue) + ylim(minValue,maxValue)

graph <- graph + theme_classic()
graph <- graph + theme(axis.text = element_text(size = 7),
                       #axis.title = element_text(size = 7, face="bold"),
                       axis.title.x = element_text(size=7, colour = "black"),
                       axis.title.y = element_text(size=7, colour = "black"),
                       text = element_text(size=7),
                       axis.text.x = element_text(angle=0,hjust=1, colour="black"),
                       axis.text.y=element_text(colour="black"),
                       axis.ticks=element_line(colour="black"),
                       #axis.title.x = element_blank(),
                       #axis.title.y = element_blank(),
                       #panel.border = element_rect(linetype = "solid", colour = "black"),
                       #panel.border = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       #panel.background = element_blank(),
                       #axis.line=element_line(colour="black"),
                       plot.title = element_text(lineheight=1.5, face="plain",size=7,hjust=0.5),
                       plot.subtitle = element_text(lineheight=1.5, face="plain",size=7,hjust=0.5),
                       legend.position="right")

print(graph)

filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/tumor_and_normal"

fileName<-paste(metabolite,"_pearson_correlation_tumor_and_normal.pdf",sep="")
fileName<-file.path(filePath,fileName)
#png(fileName,width=640,height=640,units="px")
#pdf(fileName,width=2.5,height=2.5)
pdf(fileName,width=2,height=2)
print(graph)
dev.off()


}

resultSummary<-rbind.fill(res_summary)
resultSummary<-resultSummary[order(resultSummary$pvalue),]
resultSummary$padj<-p.adjust(resultSummary$pvalue,method="BH")

filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/tumor_and_normal"

fileName<-paste("GCMS_LCMS_pearson_correlation_summary_tumor_and_normal.txt",sep="")
fileName<-file.path(filePath,fileName)

write.table(resultSummary, file=fileName, sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


######
# revise bar plot
######

res_summary<-list()

filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/tumor_and_normal/boxplot"
dir.create(filePath,recursive = TRUE)

for(metabolite in metaboliteList){
  
  #metabolite<-"pyruvate"
  
  metabolite<-"citrate"
  
  dat1<-gc_metabolite[[metabolite]]
  dat1[dat1==0]<-NA
  colnames(dat1)<-c("sampleName","metabolite_abundance")
  dat1$sampleName<-do.call(rbind,strsplit(dat1$sampleName,"_"))[,1]
  dat1$metabolite_abundance<-log2(as.numeric(dat1$metabolite_abundance))
  dat1$type<-"GC_MS"
  dat1$group<-NA
  dat1[dat1$sampleName %in% sampleName_tumor,]$group<-"tumor"
  dat1[dat1$sampleName %in% sampleName_normal,]$group<-"normal"
  
  #####
  # make relative fold change
  dat1$metabolite_abundance<-2^(dat1$metabolite_abundance)
  
  mean_tumor_vector<-mean(dat1[dat1$group %in% "tumor",]$metabolite_abundance)
  mean_normal_vector<-mean(dat1[dat1$group %in% "normal",]$metabolite_abundance)
  
  dat1$relative_fold_change<-NA
  dat1[dat1$group %in% "tumor",]$relative_fold_change<-(dat1[dat1$group %in% "tumor",]$metabolite_abundance)/(mean_normal_vector)
  dat1[dat1$group %in% "normal",]$relative_fold_change<-(dat1[dat1$group %in% "normal",]$metabolite_abundance)/(mean_normal_vector)
  
    
  #####
  
  met_vector<-metabolon_met_log2_subset[metabolite,]
  dat2<-data.frame(names(met_vector),met_vector,stringsAsFactors = FALSE)
  colnames(dat2)<-c("sampleName","metabolite_abundance")
  dat2$type<-"LC_MS"
  dat2$group<-NA
  dat2[dat2$sampleName %in% sampleName_tumor,]$group<-"tumor"
  dat2[dat2$sampleName %in% sampleName_normal,]$group<-"normal"
  
  #####
  # make relative fold change
  dat2$metabolite_abundance<-2^(dat2$metabolite_abundance)
  
  mean_tumor_vector<-mean(dat2[dat2$group %in% "tumor",]$metabolite_abundance)
  mean_normal_vector<-mean(dat2[dat2$group %in% "normal",]$metabolite_abundance)
  
  dat2$relative_fold_change<-NA
  dat2[dat2$group %in% "tumor",]$relative_fold_change<-(dat2[dat2$group %in% "tumor",]$metabolite_abundance)/(mean_normal_vector)
  dat2[dat2$group %in% "normal",]$relative_fold_change<-(dat2[dat2$group %in% "normal",]$metabolite_abundance)/(mean_normal_vector)
  
  
  
  ####
  
  plotData<-rbind(dat1,dat2)
  
  maxValue<-round(max(plotData$relative_fold_change,na.rm=TRUE))+1
  minValue<-round(min(plotData$relative_fold_change,na.rm=TRUE))-1
  
  
  #######
  #library(ggpubr)
  
  #ggscatter(plotDat,x="metabolon_abundance",y="gc_abundance",
  #          add = "reg.line", conf.int = TRUE, 
  #          cor.coef = TRUE, cor.method = "pearson",
  #          xlab = "metabolon, log2(abundance)", ylab = "GC, log2(abundance)")
  
  #######
  
  #res<-cor.test(plotDat$metabolon_abundance,plotDat$gc_abundance,method="pearson")
  #pvalue<-res$p.value
  #pvalue_for_plot<-formatC(res$p.value,format="e", digits=2)
  #coeff<-as.numeric(round(res$estimate,digit=2))
  
  #res_summary[[metabolite]]<-data.frame(metabolite,coeff,pvalue,stringsAsFactors = FALSE)
  
  
  #######
  
  #plot(log2(as.numeric(gc_pyruvate[,2])),metabolon_met_log2_subset["pyruvate",sampleName],pch=19,xlim=c(minValue,maxValue),ylim=c(minValue,maxValue))
  
  
  #plotData$group<-factor(plotData$group,levels=c("normal","tumor","tumor_class_1","tumor_class_2","tumor_class_3","tumor_class_4"))
  
  plotData$group<-factor(plotData$group,levels=c("normal","tumor"))
  plotData$type<-factor(plotData$type,levels=c("GC_MS","LC_MS"))
  
  #https://stackoverflow.com/questions/55719960/position-dodge-does-not-seem-to-work-with-stat-summary-and-x-variables-with
  
  graph <- ggplot(plotData,aes(x=type,y=relative_fold_change,fill=group))
  #graph <- graph + geom_violin()
  #graph <- graph + geom_boxplot(outlier.shape = NA)
  #graph <- graph + geom_bar(aes(fill=type),position="dodge",stat="identity")
  
  graph <- graph + stat_summary(aes(y=relative_fold_change),
                                fun="mean",geom="bar",
                                position=position_dodge(width=0.7),
                                stat="identity", 
                                width=0.6)
  graph <- graph + stat_summary(fun.data = "mean_se", 
                                geom = "errorbar",
                                position=position_dodge(width=0.7),
                                stat="identity",
                                width=0.2
                                )
  
  
  #graph <- graph + geom_jitter(position=position_jitter(0.2), size=1)
  graph <- graph + labs(title=metabolite)
  graph <- graph + xlab("Method")
  graph <- graph + ylab("Relative ratio")
  #graph <- graph + scale_color_manual(values=c("grey","red","royalblue","mediumpurple1","pink1"))
  graph <- graph + scale_fill_manual( values=c('normal' = 'grey','tumor' = 'indianred1'))
  graph <- graph + scale_y_continuous(expand = expansion(mult=c(0, 0.1),add=0))
  graph <- graph + theme_classic()
  graph <- graph + theme(axis.text = element_text(size = 7),
                         #axis.title = element_text(size = 7, face="bold"),
                         #axis.title.x = element_text(size=7),
                         axis.title.x = element_blank(),
                         axis.title.y = element_text(size=7),
                         text = element_text(size=7),
                         axis.text.x=element_text(colour="black",angle=0, hjust=0.5),
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
                         legend.text=element_text(size=7),
                         legend.key.size = unit(5,"mm"),
                         legend.position="right")
  
  print(graph)
  
  filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/tumor_and_normal/boxplot"
  
  fileName<-paste(metabolite,"_bar_plot_tumor_and_normal.pdf",sep="")
  fileName<-file.path(filePath,fileName)
  #png(fileName,width=640,height=640,units="px")
  pdf(fileName,width=2.5,height=2)
  print(graph)
  dev.off()
  
  
}


#####

res_summary<-list()

filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/tumor"
dir.create(filePath,recursive = TRUE)

for(metabolite in metaboliteList){
  
  #metabolite<-"pyruvate"
  
  dat1<-gc_metabolite[[metabolite]]
  dat1[dat1==0]<-NA
  colnames(dat1)<-c("sampleName","gc_abundance")
  
  dat1<-dat1[grep("Tumor",dat1$sampleName),]
  
  dat1$sampleName<-do.call(rbind,strsplit(dat1$sampleName,"_"))[,1]
  dat1$gc_abundance<-log2(as.numeric(dat1$gc_abundance))
  
  met_vector<-metabolon_met_log2_subset[metabolite,sampleName_tumor]
  dat2<-data.frame(names(met_vector),met_vector,stringsAsFactors = FALSE)
  colnames(dat2)<-c("sampleName","metabolon_abundance")
  
  cc<-merge(dat1,dat2,by="sampleName")
  rownames(cc)<-cc$sampleName
  cc<-cc[sampleName_tumor,]
  
  maxValue<-round(max(cc$gc_abundance,cc$metabolon_abundance,na.rm=TRUE))+1
  minValue<-round(min(cc$gc_abundance,cc$metabolon_abundance,na.rm=TRUE))-1
  
  
  plotDat<-cc
  
  #######
  #library(ggpubr)
  
  #ggscatter(plotDat,x="metabolon_abundance",y="gc_abundance",
  #          add = "reg.line", conf.int = TRUE, 
  #          cor.coef = TRUE, cor.method = "pearson",
  #          xlab = "metabolon, log2(abundance)", ylab = "GC, log2(abundance)")
  
  #######
  
  res<-cor.test(plotDat$metabolon_abundance,plotDat$gc_abundance,method="pearson")
  pvalue<-res$p.value
  pvalue_for_plot<-formatC(res$p.value,format="e", digits=2)
  coeff<-as.numeric(round(res$estimate,digit=2))
  
  res_summary[[metabolite]]<-data.frame(metabolite,coeff,pvalue,stringsAsFactors = FALSE)
  
  
  #######
  
  #plot(log2(as.numeric(gc_pyruvate[,2])),metabolon_met_log2_subset["pyruvate",sampleName],pch=19,xlim=c(minValue,maxValue),ylim=c(minValue,maxValue))
  
  
  graph <- ggplot(plotDat,aes(x=metabolon_abundance,y=gc_abundance))
  graph <- graph + geom_point(size=1,shape=19)
  graph <- graph + labs(title=metabolite)
  graph <- graph + labs(subtitle=paste("R = ",coeff,", ","p = ",pvalue_for_plot,sep=""))
  graph <- graph + geom_smooth(method=lm)
  graph <- graph + xlab("log2(LC-MS abundance)")
  graph <- graph + ylab("log2(GC-MS abundance)")
  
  #graph <- graph + geom_abline(intercept = 0,slope=1,color="grey")
  #graph <- graph + xlim(minValue,maxValue) + ylim(minValue,maxValue)
  
  graph <- graph + theme_classic()
  graph <- graph + theme(axis.text = element_text(size = 7),
                         #axis.title = element_text(size = 7, face="bold"),
                         axis.title.x = element_text(size=7, colour = "black"),
                         axis.title.y = element_text(size=7, colour = "black"),
                         text = element_text(size=7),
                         axis.text.x = element_text(angle=0,hjust=1, colour="black"),
                         axis.text.y=element_text(colour="black"),
                         axis.ticks=element_line(colour="black"),
                         #axis.title.x = element_blank(),
                         #axis.title.y = element_blank(),
                         #panel.border = element_rect(linetype = "solid", colour = "black"),
                         #panel.border = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         #panel.background = element_blank(),
                         #axis.line=element_line(colour="black"),
                         plot.title = element_text(lineheight=1.5, face="plain",size=7,hjust=0.5),
                         plot.subtitle = element_text(lineheight=1.5, face="plain",size=7,hjust=0.5),
                         legend.position="right")
  
  print(graph)
  
  filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/tumor"
  
  fileName<-paste(metabolite,"_pearson_correlation_tumor.pdf",sep="")
  fileName<-file.path(filePath,fileName)
  #png(fileName,width=640,height=640,units="px")
  pdf(fileName,width=2.5,height=2.5)
  print(graph)
  dev.off()
  
  
}

resultSummary<-rbind.fill(res_summary)
resultSummary<-resultSummary[order(resultSummary$pvalue),]
resultSummary$padj<-p.adjust(resultSummary$pvalue,method="BH")

filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/tumor"

fileName<-paste("GCMS_LCMS_pearson_correlation_summary_tumor.txt",sep="")
fileName<-file.path(filePath,fileName)

write.table(resultSummary, file=fileName, sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

#####

res_summary<-list()

filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/normal"
dir.create(filePath,recursive = TRUE)

for(metabolite in metaboliteList){
  
  #metabolite<-"pyruvate"
  
  dat1<-gc_metabolite[[metabolite]]
  dat1[dat1==0]<-NA
  colnames(dat1)<-c("sampleName","gc_abundance")
  
  dat1<-dat1[grep("Normal",dat1$sampleName),]
  
  dat1$sampleName<-do.call(rbind,strsplit(dat1$sampleName,"_"))[,1]
  dat1$gc_abundance<-log2(as.numeric(dat1$gc_abundance))
  
  met_vector<-metabolon_met_log2_subset[metabolite,sampleName_normal]
  dat2<-data.frame(names(met_vector),met_vector,stringsAsFactors = FALSE)
  colnames(dat2)<-c("sampleName","metabolon_abundance")
  
  cc<-merge(dat1,dat2,by="sampleName")
  rownames(cc)<-cc$sampleName
  cc<-cc[sampleName_normal,]
  
  maxValue<-round(max(cc$gc_abundance,cc$metabolon_abundance,na.rm=TRUE))+1
  minValue<-round(min(cc$gc_abundance,cc$metabolon_abundance,na.rm=TRUE))-1
  
  
  plotDat<-cc
  
  #######
  #library(ggpubr)
  
  #ggscatter(plotDat,x="metabolon_abundance",y="gc_abundance",
  #          add = "reg.line", conf.int = TRUE, 
  #          cor.coef = TRUE, cor.method = "pearson",
  #          xlab = "metabolon, log2(abundance)", ylab = "GC, log2(abundance)")
  
  #######
  
  res<-cor.test(plotDat$metabolon_abundance,plotDat$gc_abundance,method="pearson")
  pvalue<-res$p.value
  pvalue_for_plot<-formatC(res$p.value,format="e", digits=2)
  coeff<-as.numeric(round(res$estimate,digit=2))
  
  res_summary[[metabolite]]<-data.frame(metabolite,coeff,pvalue,stringsAsFactors = FALSE)
  
  
  #######
  
  #plot(log2(as.numeric(gc_pyruvate[,2])),metabolon_met_log2_subset["pyruvate",sampleName],pch=19,xlim=c(minValue,maxValue),ylim=c(minValue,maxValue))
  
  
  graph <- ggplot(plotDat,aes(x=metabolon_abundance,y=gc_abundance))
  graph <- graph + geom_point(size=1,shape=19)
  graph <- graph + labs(title=metabolite)
  graph <- graph + labs(subtitle=paste("R = ",coeff,", ","p = ",pvalue_for_plot,sep=""))
  graph <- graph + geom_smooth(method=lm)
  graph <- graph + xlab("log2(LC-MS abundance)")
  graph <- graph + ylab("log2(GC-MS abundance)")
  
  #graph <- graph + geom_abline(intercept = 0,slope=1,color="grey")
  #graph <- graph + xlim(minValue,maxValue) + ylim(minValue,maxValue)
  
  graph <- graph + theme_classic()
  graph <- graph + theme(axis.text = element_text(size = 7),
                         #axis.title = element_text(size = 7, face="bold"),
                         axis.title.x = element_text(size=7, colour = "black"),
                         axis.title.y = element_text(size=7, colour = "black"),
                         text = element_text(size=7),
                         axis.text.x = element_text(angle=0,hjust=1, colour="black"),
                         axis.text.y=element_text(colour="black"),
                         axis.ticks=element_line(colour="black"),
                         #axis.title.x = element_blank(),
                         #axis.title.y = element_blank(),
                         #panel.border = element_rect(linetype = "solid", colour = "black"),
                         #panel.border = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         #panel.background = element_blank(),
                         #axis.line=element_line(colour="black"),
                         plot.title = element_text(lineheight=1.5, face="plain",size=7,hjust=0.5),
                         plot.subtitle = element_text(lineheight=1.5, face="plain",size=7,hjust=0.5),
                         legend.position="right")
  
  print(graph)
  
  filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/normal"
  
  fileName<-paste(metabolite,"_pearson_correlation_normal.pdf",sep="")
  fileName<-file.path(filePath,fileName)
  #png(fileName,width=640,height=640,units="px")
  pdf(fileName,width=2.5,height=2.5)
  print(graph)
  dev.off()
  
  
}

resultSummary<-rbind.fill(res_summary)
resultSummary<-resultSummary[order(resultSummary$pvalue),]
resultSummary$padj<-p.adjust(resultSummary$pvalue,method="BH")

filePath<-"~/work/Ed_lab/HCC_project/results/associate_GC_with_Metabolon/06_23_2021/normal"

fileName<-paste("GCMS_LCMS_pearson_correlation_summary_normal.txt",sep="")
fileName<-file.path(filePath,fileName)

write.table(resultSummary, file=fileName, sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

#####




