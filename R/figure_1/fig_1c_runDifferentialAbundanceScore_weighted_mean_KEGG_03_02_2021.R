library(plyr)
library(data.table)

######
# note: add fold change as weight to differential abundance score
######

date<-"03_02_2021"


#######
# handle differential metabolite abundance list
#######

filePath<-"/Users/mil20411/work/Ed_lab/HCC_project/data/metabotools_HCC"
fileName<-"HCC_metabolite_info_04_22_2020.Rda"
fileName<-file.path(filePath,fileName)

load(file=fileName)

metinfoSubset<-metinfo[,c("name","H_KEGG")]


######

testCaseNameList<-c("tumor_vs_normal",
                    "hcc_tumor_vs_hcc_normal",
                    "hcc_vs_pd",
                    "hcc_vs_ptctv",
                    #"hcc_HWIDE_vs_hcc_HMIN", # there is no significantly different metabolites in this situation 
                    "hcc_vs_pd_and_ptctv",
                    "hcc_vs_ha",
                    "pd_ptctv_tumor_vs_pd_ptctv_normal",
                    "ha_pd_ptctv_tumor_vs_ha_pd_ptctv_normal")

#testCaseName<-"hcc_vs_pd_and_ptctv"

for(j in 1:length(testCaseNameList)){

    testCaseName<-testCaseNameList[j]
    message(sprintf("Work on %s",testCaseName))

    filePath<-file.path("/Users/mil20411/work/Ed_lab/HCC_project/results/diffMetaboliteAbundance/02_09_2021")
    #fileName<-"HCC_differential_metabolite_abundance_wilcox_test_tumor_vs_normal.txt"
    #fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_vs_pd.txt"
    #fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_vs_ptctv.txt"
    #fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_HWIDE_vs_hcc_HMIN.txt"
    #fileName<-"HCC_differential_metabolite_abundance_wilcox_test_hcc_vs_pd_and_ptctv.txt"
    
    fileName<-paste("HCC_differential_metabolite_abundance_wilcox_test_",testCaseName,".txt",sep="")
    
    fileName<-file.path(filePath,fileName)
    
    # Be careful about line reading lost
    dat<-fread(file=fileName,sep="\t",header=TRUE,data.table = FALSE)
    #dat<-read.table(file=fileName,sep="\t",header=TRUE,stringsAsFactors = FALSE,fill=" ")
    #dat$name<-rownames(dat)
    
    #matchedDiffMetAbundance<-merge(metinfoSubset,dat,by="name",all.x=TRUE)
    matchedDiffMetAbundance<-merge(metinfo,dat,by="name",all.x=TRUE)
    
    matchedDiffMetAbundance$rankValue<-matchedDiffMetAbundance$log2FC*(-log10(matchedDiffMetAbundance$padj))
    
    matchedDiffMetAbundanceResult<-matchedDiffMetAbundance[!is.na(matchedDiffMetAbundance$H_KEGG) & !is.na(matchedDiffMetAbundance$padj),]
    #matchedDiffMetAbundance<-matchedDiffMetAbundance[!is.na(matchedDiffMetAbundance$H_KEGG) & !is.na(matchedDiffMetAbundance$rankValue),]
    
    
    # remove duplicated entries by checking KEGG ID 
    matchedDiffMetAbundanceResult<-matchedDiffMetAbundanceResult[!duplicated(matchedDiffMetAbundanceResult$H_KEGG),]
    
    #threshold<-0.05
    #matchedDiffMetAbundance<-matchedDiffMetAbundanceResult[matchedDiffMetAbundanceResult$padj<threshold,]
    
    
    matchedDiffMetAbundanceList<-matchedDiffMetAbundance$log2FC
    names(matchedDiffMetAbundanceList)<-matchedDiffMetAbundance$H_KEGG
    
    #matchedMetAbundanceRank<- -log10(matchedDiffMetAbundance$padj)
    matchedMetAbundanceRank<-matchedDiffMetAbundance$rankValue
    names(matchedMetAbundanceRank)<-matchedDiffMetAbundance$H_KEGG
    
    #matchedMetAbundanceRank<-sort(matchedMetAbundanceRank,decreasing = FALSE)
    
    matchedMetAbundanceRank<-sort(matchedMetAbundanceRank,decreasing = TRUE)
    
    #######
    # calculate differential abundance score on differential metabolite abundance
    #######
    
    filePath<-"~/work/Ed_lab/HCC_project/data/KEGG_metabolic_pathway/parsed_data"
    fileName<-paste("KEGG_metabolic_pathway_metaboliteInfoAll.txt",sep="")
    fileName<-file.path(filePath,fileName)
    
    pathwayDat<-fread(file=fileName,header=TRUE,sep="\t",stringsAsFactors = FALSE,data.table=FALSE)
    
    
    pathwayList<-unique(pathwayDat$pathwayName)
    
    DA_results<-list()
    
    for(idx in 1:length(pathwayList)){
      
      pathwayName<-pathwayList[idx]
      
      message(sprintf("work on %s %s/%s",pathwayName,idx,length(pathwayList)))
      
      metaboliteList<-pathwayDat[pathwayDat$pathwayName %in% pathwayName,]$keggMetaboliteID
      
      metaboliteListMeasured<-matchedDiffMetAbundanceResult[matchedDiffMetAbundanceResult$H_KEGG %in% metaboliteList,]$H_KEGG
      
      numOfMetabolites<-length(metaboliteListMeasured)
      
      ##### old way to calculate differential abundance score
      if(FALSE){
        threshold<-0.05
        sig_metabolite_in_pathway<-matchedDiffMetAbundanceResult[(matchedDiffMetAbundanceResult$H_KEGG %in% metaboliteListMeasured) & matchedDiffMetAbundanceResult$padj<threshold,]
        
        numOfMetaboliteUp<-sum(sig_metabolite_in_pathway$log2FC>0)
        numOfMetaboliteDown<-sum(sig_metabolite_in_pathway$log2FC<0)
        
        differentialAbundanceScore<-(numOfMetaboliteUp-numOfMetaboliteDown)/numOfMetabolites
      }
      
      ##### 11/09/2020 revision to calculate differential abundance score
      if(FALSE){
        threshold<-0.05
        sig_metabolite_in_pathway<-matchedDiffMetAbundanceResult[(matchedDiffMetAbundanceResult$H_KEGG %in% metaboliteListMeasured) & matchedDiffMetAbundanceResult$padj<threshold,]
        
        numOfMetaboliteUp<-sum(sig_metabolite_in_pathway$log2FC>0)
        numOfMetaboliteDown<-sum(sig_metabolite_in_pathway$log2FC<0)
        
        DAscoreByFoldChange<-sum(sig_metabolite_in_pathway$log2FC)
        
        differentialAbundanceScore<-(DAscoreByFoldChange)/numOfMetabolites
      }
      #####
      ##### 12/22/2020 revision 
      ##### use weighted mean to calculate differential abundance score
      
      if(TRUE){
        
        
        threshold<-0.05
        
        sig_metabolite_in_pathway<-matchedDiffMetAbundanceResult[(matchedDiffMetAbundanceResult$H_KEGG %in% metaboliteListMeasured) & matchedDiffMetAbundanceResult$padj<threshold,]
        
        #measured_metabolite_in_pathway<-matchedDiffMetAbundanceResult[(matchedDiffMetAbundanceResult$H_KEGG %in% metaboliteListMeasured),]
                                                                 
        numOfMetaboliteUp<-sum(sig_metabolite_in_pathway$log2FC>0)
        numOfMetaboliteDown<-sum(sig_metabolite_in_pathway$log2FC<0)
        
        #pathwayFrame<-matchedDiffMetAbundanceResult[(matchedDiffMetAbundanceResult$H_KEGG %in% metaboliteListMeasured),]
        pathwayFrame<-matchedDiffMetAbundanceResult[(matchedDiffMetAbundanceResult$H_KEGG %in% metaboliteListMeasured),]
                      
        if(numOfMetabolites>0){
          
          pathwayFrame$explainVector<-0
          
          if( sum(pathwayFrame$padj<threshold & pathwayFrame$log2FC>0) > 0 ){
             pathwayFrame[pathwayFrame$padj<threshold & pathwayFrame$log2FC>0,]$explainVector<-1
          }
          
          if( sum(pathwayFrame$padj<threshold & pathwayFrame$log2FC<0) > 0 ){
             pathwayFrame[pathwayFrame$padj<threshold & pathwayFrame$log2FC<0,]$explainVector<-(-1)
          }
            explainVector<-pathwayFrame$explainVector
            weightVector<-abs(pathwayFrame$log2FC)/sum(abs(pathwayFrame$log2FC))
            #weightVector<-abs(pathwayFrame$log2FC)
            #euclideanNorm<-sqrt(sum(weightVector^2))
            #weightVectorNormalized<-weightVector/euclideanNorm
            
            #differentialAbundanceScore<-sum(explainVector * weightVectorNormalized
            
            # weighted mean 
            # https://stat.ethz.ch/R-manual/R-devel/library/stats/html/weighted.mean.html
            differentialAbundanceScore<-weighted.mean(explainVector,weightVector)
            
            
        }else{
          differentialAbundanceScore<-0
        }
      
      }
      #####
      
      DA_results[[pathwayName]]<-data.frame(pathwayName,numOfMetaboliteUp,numOfMetaboliteDown,numOfMetabolites,differentialAbundanceScore,stringsAsFactors = FALSE)
      
    }
    
    DA_results<-rbind.fill(DA_results)
    DA_results<-DA_results[order(-DA_results$differentialAbundanceScore),]
    
    filePath<-file.path("/Users/mil20411/work/Ed_lab/HCC_project/results/pathway_enrichment/differential_abundance_score",date,"weighted_mean")
    dir.create(filePath, showWarnings = TRUE, recursive = TRUE)
    fileName<-paste("HCC_metabolite_differential_abundance_score_weighted_mean_KEGG_pathway_",testCaseName,"_",date,".txt",sep="")
    
    fileName<-file.path(filePath,fileName)
    
    write.table(DA_results,file=fileName,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)
    
    ########
    
    library(ggplot2)
    
    DA_results<-DA_results[(DA_results$numOfMetabolites !=0),]
    DA_results$pathwaySize<-DA_results$numOfMetabolites
    
    DA_results$status<-"none"
    
    if(sum(DA_results$differentialAbundanceScore>0)>0){
    DA_results[DA_results$differentialAbundanceScore>0,]$status<-"up"
    }
    
    if(sum(DA_results$differentialAbundanceScore<0)>0){
    DA_results[DA_results$differentialAbundanceScore<0,]$status<-"down"
    }
    
    
    ### order differential abundance score
    
    criteria1 <- DA_results$differentialAbundanceScore>0
    subsetDatUp<-DA_results[criteria1,]
    
    criteria2 <- order(subsetDatUp$differentialAbundanceScore,subsetDatUp$pathwaySize)
    
    descOrderUp<-subsetDatUp[ criteria2, ]$pathwayName
    
    subsetDatDown<-DA_results[!criteria1,]
    criteria3 <- order(subsetDatDown$differentialAbundanceScore,-subsetDatDown$pathwaySize) 
    
    descOrderDown<-subsetDatDown[ criteria3, ]$pathwayName
    
    descOrder<-c(descOrderDown, descOrderUp)
    
    #DA_results$differentialAbundanceScore<-factor(DA_results$differentialAbundanceScore)
    
    DA_results$pathwayName<-factor(DA_results$pathwayName, 
                                   levels=descOrder)
    
    #DA_results<-DA_results[order(-DA_results$differentialAbundanceScore),]
    
    ####
    
    # only show pathway with more than 5 measured metabolites
    pathwaySize_threshold<-5
    
    DA_results_removed<-DA_results[DA_results$pathwaySize<pathwaySize_threshold,]
    DA_results_kept<-DA_results[DA_results$pathwaySize>=pathwaySize_threshold,]
    
    
    ####
    
    #myColor<-c("none"="grey","up"="red","down"="blue")
    mid<-0
    
    DA_score_range<-ceiling(max(max(DA_results_kept$differentialAbundanceScore),abs(min(DA_results_kept$differentialAbundanceScore))))
    
    #graph <- ggplot(data=DA_results, aes(x=differentialAbundanceScore, y=pathwayName, color=differentialAbundanceScore))
    graph <- ggplot(data=DA_results_kept, aes(x=differentialAbundanceScore, y=pathwayName, color=differentialAbundanceScore))
    #graph <- ggplot(data=DA_results_removed, aes(x=differentialAbundanceScore, y=pathwayName, color=differentialAbundanceScore))
    graph <- graph + geom_segment( aes(x=0, xend=differentialAbundanceScore,y=pathwayName, yend=pathwayName,color=differentialAbundanceScore ))
    graph <- graph + geom_point( aes(size=pathwaySize), fill=alpha("orange", 0.3), alpha=1)
    graph <- graph + scale_size(range=c(0.5,3))
    #graph <- graph + expand_limits(x=c(-1,1))
    graph <- graph + expand_limits(x=c(-DA_score_range,DA_score_range))
    #graph <- graph + scale_color_manual(values=myColor)
    graph <- graph + scale_color_gradient2(midpoint=mid, 
                                           low="blue",
                                           mid="grey",
                                           high="red",
                                           space="Lab")
    graph <- graph + geom_vline(xintercept = 0, linetype="dotted",col="grey50")
    #graph <- graph  + ylab("Pathways") 
    graph <- graph  + ylab("") 
    graph <- graph + xlab("Differential Abundance Score")
    graph <- graph + scale_y_discrete(expand = expansion(add=1))
    #graph <- graph + coord_flip()
    graph <- graph + theme_classic()
    graph <- graph + theme(  
                             #axis.title = element_text(size = 7, face="bold"),
                             axis.title.x = element_text(size=7),
                             axis.title.y = element_text(size=7),
                             #text = element_text(size=12),
                             #axis.text.x = element_text(), 
                             #axis.title.x = element_blank(),
                             #axis.title.y = element_blank(),
                             #panel.border = element_rect(linetype = "solid", colour = "black"),
                             #panel.border = element_blank(),
                             #panel.grid.major = element_blank(), 
                             #panel.grid.minor = element_blank(),
                             panel.grid.major.y = element_line(colour = "grey90"),
                             panel.grid.minor.y = element_line(colour = "grey90"),
                             panel.grid.major.x = element_line(colour = "grey90"),
                             panel.grid.minor.x =  element_blank(),
                             axis.text.x=element_text(colour="black", size = 7, angle= 0, hjust = 1),
                             axis.text.y=element_text(colour="black", size = 7),
                             axis.ticks=element_line(colour="black"),
                             axis.line=element_line(colour="black"),
                             #panel.background = element_blank(),
                             #axis.line=element_line(colour="black",size=0.7),
                             plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                             plot.subtitle = element_text(lineheight=1.5, face="plain",size=7,hjust=0.5),
                             plot.margin= margin(10, 10, 3, 1, "pt"),
                             legend.position="none")
    
    print(graph)
    
    filePath<-file.path("/Users/mil20411/work/Ed_lab/HCC_project/results/pathway_enrichment/differential_abundance_score",date,"weighted_mean")
    dir.create(filePath, showWarnings = TRUE, recursive = TRUE)
    #fileName<-"HCC_differential_abundance_score.pdf"
    fileName<-paste("HCC_meatbolite_DA_score_weighted_mean_KEGG_pathway_removed_less_than_5_measured_metabolites_in_pathway_",testCaseName,"_",date,".pdf",sep="")
    #fileName<-"HCC_differential_abundance_score_show_less_than_5_measured_metabolites_in_pathway.pdf"
    
    fileName<-file.path(filePath,fileName)
    
    pdf(fileName,height=5.2,width=3.5)
    print(graph)
    dev.off()

#######
    
    graph <- ggplot(data=DA_results_kept, aes(x=differentialAbundanceScore, y=pathwayName, color=differentialAbundanceScore))
    #graph <- ggplot(data=DA_results_removed, aes(x=differentialAbundanceScore, y=pathwayName, color=differentialAbundanceScore))
    graph <- graph + geom_segment( aes(x=0, xend=differentialAbundanceScore,y=pathwayName, yend=pathwayName,color=differentialAbundanceScore ))
    graph <- graph + geom_point( aes(size=pathwaySize), fill=alpha("orange", 0.3), alpha=1)
    graph <- graph + scale_size(range=c(0.5,3))
    #graph <- graph + expand_limits(x=c(-1,1))
    graph <- graph + expand_limits(x=c(-DA_score_range,DA_score_range))
    #graph <- graph + scale_color_manual(values=myColor)
    graph <- graph + scale_color_gradient2(midpoint=mid, 
                                           low="blue",
                                           mid="grey",
                                           high="red",
                                           space="Lab")
    graph <- graph + geom_vline(xintercept = 0, linetype="dotted",col="grey50")
    graph <- graph  + ylab("Pathways") 
    graph <- graph + xlab("Differential Abundance Score")
    graph <- graph + scale_y_discrete(expand = expansion(add=1))
    #graph <- graph + coord_flip()
    graph <- graph + theme_classic()
    graph <- graph + theme(  
      #axis.title = element_text(size = 7, face="bold"),
      axis.title.x = element_text(size=7),
      axis.title.y = element_text(size=7),
      #text = element_text(size=12),
      #axis.text.x = element_text(), 
      #axis.title.x = element_blank(),
      #axis.title.y = element_blank(),
      #panel.border = element_rect(linetype = "solid", colour = "black"),
      #panel.border = element_blank(),
      #panel.grid.major = element_blank(), 
      #panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      panel.grid.minor.y = element_line(colour = "grey90"),
      panel.grid.major.x = element_line(colour = "grey90"),
      panel.grid.minor.x =  element_blank(),
      axis.text.x=element_text(colour="black", size = 7, angle= 0, hjust = 1),
      axis.text.y=element_text(colour="black", size = 7),
      axis.ticks=element_line(colour="black"),
      axis.line=element_line(colour="black"),
      #panel.background = element_blank(),
      #axis.line=element_line(colour="black",size=0.7),
      plot.title = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
      plot.subtitle = element_text(lineheight=1.5, face="plain",size=7,hjust=0.5),
      plot.margin= margin(10, 10, 3, 3, "pt"),
      legend.position="right")
    
    print(graph)
    
    filePath<-file.path("/Users/mil20411/work/Ed_lab/HCC_project/results/pathway_enrichment/differential_abundance_score",date,"weighted_mean")
    dir.create(filePath, showWarnings = TRUE, recursive = TRUE)
    #fileName<-"HCC_differential_abundance_score.pdf"
    fileName<-paste("HCC_meatbolite_DA_score_weighted_mean_KEGG_pathway_removed_less_than_5_measured_metabolites_in_pathway_",testCaseName,"_",date,"_legend.pdf",sep="")
    #fileName<-"HCC_differential_abundance_score_show_less_than_5_measured_metabolites_in_pathway.pdf"
    
    fileName<-file.path(filePath,fileName)
    
    pdf(fileName,height=5.2,width=6)
    print(graph)
    dev.off()
    
    
}







