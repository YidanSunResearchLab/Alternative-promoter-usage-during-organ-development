prewd="Mammalian2019Promoter/CM2019RNAseqMouse"
stage.order = c("s10_Female","s10_Male","s11_Female","s11_Male","s12_Female","s12_Male","s13_Female","s13_Male","s14_Female","s14_Male",
                "s15_Female","s15_Male","s16_Female","s16_Male","s17_Female","s17_Male","s18_Female","s18_Male","s0dpb_Female","s0dpb_Male",
                "s3dpb_Female","s3dpb_Male","s2wpb_Female","s2wpb_Male","s4wpb_Female","s4wpb_Male","s9wpb_Female","s9wpb_Male")
sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects","sample.sheet.txt"),header=TRUE, row.names = 1, stringsAsFactors = FALSE)


library(GenomicRanges)
library(GenomicFeatures)
library(ensembldb)
library(GenomicAlignments)
library(AnnotationDbi)
library(dplyr)
library(proActiv)
library(ggplot2)
library(reshape2)
# library(preprocessCore) #function normalize.quantiles
library(Rtsne)
library(Seurat)
library(pheatmap)
library(doParallel)
library(foreach)

method="Proactiv"
# method="Dexseq"
# norm.method = "edger"
norm.method = "deseq2"
gender = "Male"
# gender = "Female"
# gender = "ale"

tissue.order = c("Brain", "Cerebellum", "Heart", "Kidney" , "Liver", "Ovary", "Testis")
mc.cores= 20
basic.objects=c("basic.objects",ls())

message(method)
message(norm.method)
################################################
###TSNE plot for all tissues and development stages
################################################
if(TRUE){
  dir.create(file.path(prewd,"Output","TissueSpecific"))
  dir.create(file.path(prewd,"Output","TissueSpecific",method))
  dir.create(file.path(prewd,"Output","TissueSpecific",method,norm.method))
  
  ############################################################################
  ##calculate standard variation and use the top 2000 most variable promoters
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.RData'))
  absolute.promoter.activity.star.mean.sd = apply(absolute.promoter.activity.star.mean,1,sd)
  absolute.promoter.activity.star.mean.sd = sort(absolute.promoter.activity.star.mean.sd,decreasing = TRUE)[1:2000]
  
  ##TSNE calculation
  absolute.promoter.activity.star.mean.selected = absolute.promoter.activity.star.mean[names(absolute.promoter.activity.star.mean.sd),]
  set.seed(100) 
  absolute.promoter.activity.star.mean.top2000.tsne <- Rtsne(t(absolute.promoter.activity.star.mean.selected),check_duplicates=FALSE,theta=0.5,num_threads=mc.cores) # Run TSNE
  save(absolute.promoter.activity.star.mean.top2000.tsne,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.top2000.tsne.RData'))
  
  # library(umap)
  # set.seed(100) # Sets seed for reproducibility
  # umap.out.basic <- umap::umap(gene.states.dataframe) # Run TSNE
  # save(umap.out.basic,file=file.path("TTseqObjects",directory, paste0(methodname,".umap.out.basic.",gsub("[*.]","",sample),".",binned,".RData")))
  
  ############################################################################
  ##calculate standard variation and use the top 2000 most variable genes
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'gene.expression.star.mean.RData'))
  gene.expression.star.mean.sd = apply(gene.expression.star.mean,1,sd)
  gene.expression.star.mean.sd = sort(gene.expression.star.mean.sd,decreasing = TRUE)[1:2000]
  
  ##TSNE calculation
  gene.expression.star.mean.selected = gene.expression.star.mean[names(gene.expression.star.mean.sd),]
  set.seed(110) 
  gene.expression.star.mean.top2000.tsne <- Rtsne(t(gene.expression.star.mean.selected),check_duplicates=FALSE,theta=0.5,num_threads=mc.cores) # Run TSNE
  save(gene.expression.star.mean.top2000.tsne,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'gene.expression.star.mean.top2000.tsne.RData'))
  
}

##Plot TSNE
if(TRUE){
  build.tsne.plot = function(tsne.out,genes,name){
    tsne.plot.data.frame = as.data.frame(tsne.out$Y) #
    colnames(tsne.plot.data.frame) = c("tSNE1","tSNE2")
    rownames(tsne.plot.data.frame) = colnames(genes)
    tsne.plot.data.frame$tissue = gsub("_.*","",rownames(tsne.plot.data.frame))
    tsne.plot.data.frame$stage = gsub(".*_P.*","LateStage",gsub(".*_E.*","EarlyStage",gsub("pmonth","LateStage",rownames(tsne.plot.data.frame)))) #gsub("[0-9]*","",gsub(".*_","", gsub(".*_","",rownames(tsne.plot.data.frame)) #
    tsne.plot.data.frame$sex = gsub(".*_","",rownames(tsne.plot.data.frame)) 
    if(length(grep("Both",name))==0){
      tsne.plot.data.frame = tsne.plot.data.frame[tsne.plot.data.frame$sex==gsub(".*_", "",name),]
    } else {
      tsne.plot.data.frame$stage = tsne.plot.data.frame$sex
    }
    
    sample.color = c("#323edd","#364f6b","#fe346e","#f2a51a","#0c9463","#b61aae","#ea6227")
    p = ggplot(tsne.plot.data.frame, aes(tSNE1,tSNE2,color=tissue,fill=tissue, shape=stage)) +
      # stat_density2d(aes(fill=state, alpha = ..level..),geom="polygon") +
      geom_point(alpha = 0.9, size = 2) +
      scale_color_manual(values = sample.color) +
      scale_fill_manual(values = sample.color) +
      # scale_shape_manual(values=c(1,2)) +
      # ggtitle(paste0("t-SNE plot")) +
      theme_bw(base_size = 14) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(color="black"))
    print(p)
    ggsave(p,file=file.path(prewd,"Output","TissueSpecific",method,norm.method,paste0('all.samples.',name,'.tsne.pdf')),width=6,height = 4,units = "in",dpi=300)
  }
  
  ## promoters
  ##calculate standard variation and use the top 2000 most variable promoters
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.top2000.tsne.RData'))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.RData'))
  absolute.promoter.activity.star.mean.sd = apply(absolute.promoter.activity.star.mean,1,sd)
  absolute.promoter.activity.star.mean.sd = sort(absolute.promoter.activity.star.mean.sd,decreasing = TRUE)[1:2000]
  absolute.promoter.activity.star.mean.selected = absolute.promoter.activity.star.mean[names(absolute.promoter.activity.star.mean.sd),]
  
  build.tsne.plot(tsne.out=absolute.promoter.activity.star.mean.top2000.tsne, genes=absolute.promoter.activity.star.mean.selected, name="promoter2000_Both") #Fig. 1g
  build.tsne.plot(tsne.out=absolute.promoter.activity.star.mean.top2000.tsne, genes=absolute.promoter.activity.star.mean.selected, name="promoter2000_Male") #Fig. 1h
  build.tsne.plot(tsne.out=absolute.promoter.activity.star.mean.top2000.tsne, genes=absolute.promoter.activity.star.mean.selected, name="promoter2000_Female") #Fig. 1i
  
  ## genes
  ##calculate standard variation and use the top 2000 most variable promoters
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'gene.expression.star.mean.top2000.tsne.RData'))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'gene.expression.star.mean.RData'))
  gene.expression.star.mean.sd = apply(gene.expression.star.mean,1,sd)
  gene.expression.star.mean.sd = sort(gene.expression.star.mean.sd,decreasing = TRUE)[1:2000]
  gene.expression.star.mean.selected = gene.expression.star.mean[names(gene.expression.star.mean.sd),]
  
  build.tsne.plot(tsne.out=gene.expression.star.mean.top2000.tsne, genes=gene.expression.star.mean.selected, name="gene2000_Both") #Fig. S1b
  build.tsne.plot(tsne.out=gene.expression.star.mean.top2000.tsne, genes=gene.expression.star.mean.selected, name="gene2000_Male") #Fig. S1c
  build.tsne.plot(tsne.out=gene.expression.star.mean.top2000.tsne, genes=gene.expression.star.mean.selected, name="gene2000_Female") #Fig. S1d
  
}

################################################
###Tissue specifc analysis with Seurat
################################################

if(TRUE){
  dir.create(file.path(prewd, "ProcessedData", "TissueObjects"))
  dir.create(file.path(prewd, "ProcessedData", "TissueObjects",method))
  dir.create(file.path(prewd, "ProcessedData", "TissueObjects",method,norm.method))
  dir.create(file.path(prewd, "Output", "TissueSpecific"))
  dir.create(file.path(prewd, "Output", "TissueSpecific",method))
  dir.create(file.path(prewd, "Output", "TissueSpecific",method,norm.method))
  
  #load counts
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.RData'))
  colnames(absolute.promoter.activity.star.mean) = gsub("_","-",colnames(absolute.promoter.activity.star.mean))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivity.star.RData"))
  colnames(absolutePromoterActivity.star) = gsub("_","-",colnames(absolutePromoterActivity.star))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"normalizedPromoterCounts.star.RData"))
  colnames(normalizedPromoterCounts.star) = gsub("_","-",colnames(normalizedPromoterCounts.star))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,'promoterCounts.star.RData'))
  colnames(promoterCounts.star) = gsub("_","-",colnames(promoterCounts.star))
  dcounts = absolutePromoterActivity.star[,-1:-2] #normalizedPromoterCounts.star #absolute.promoter.activity.star.mean #
  tissue.raw.counts.object = CreateSeuratObject(counts = dcounts, project = "HumanTissue", min.cells = 1, min.features = 0)
  pdf(file.path(prewd, "Output", "TissueSpecific",method,norm.method, "seurat.feature.vlnplot.pdf" ))
  print(VlnPlot(tissue.raw.counts.object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) )
  dev.off()
  
  ##########################
  ##Normalization
  ##########################
  # tissue.raw.counts.object <- NormalizeData(tissue.raw.counts.object, normalization.method = "LogNormalize", scale.factor = 10000)

  ##########################
  ##feature selection (highly variable features)
  ##########################
  tissue.raw.counts.object <- FindVariableFeatures(tissue.raw.counts.object, selection.method = "vst", nfeatures = 2000)
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(tissue.raw.counts.object)
  plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(tissue.raw.counts.object), 20), repel = TRUE)
  pdf(file.path(prewd, "Output", "TissueSpecific",method,norm.method, "seurat.feature.selection.pdf" ))
  print(plot2)
  dev.off()
  
  ##########################################################################
  ##Scaling the data and removing unwanted sources of variation
  ##########################################################################
  tissue.raw.counts.object <- ScaleData(object = tissue.raw.counts.object, vars.to.regress = c("nCount_RNA")) #features = rownames(tissue.raw.counts.object), 
  
  ########################################
  ##Perform linear dimensional reduction
  ########################################
  tissue.raw.counts.object <- RunPCA(tissue.raw.counts.object, features = VariableFeatures(object = tissue.raw.counts.object)) #, npcs = 5
  print(tissue.raw.counts.object[["pca"]], dims = 1:5, nfeatures = 5)

  pdf(file.path(prewd, "Output", "TissueSpecific", method,norm.method, "seurat.PCA.runPCA.pdf" ))
  print(VizDimLoadings(tissue.raw.counts.object, dims = 1:2, reduction = "pca"))
  print(DimPlot(tissue.raw.counts.object, reduction = "pca"))
  print(DimHeatmap(tissue.raw.counts.object, dims = 1:5, cells = ncol(dcounts), balanced = TRUE))
  dev.off()

  ########################################
  ##Determine the ‘dimensionality’ of the dataset (PC number)
  ########################################
  tissue.raw.counts.object <- JackStraw(tissue.raw.counts.object, num.replicate = 100)
  tissue.raw.counts.object <- ScoreJackStraw(tissue.raw.counts.object, dims = 1:15)
  pdf(file.path(prewd, "Output", "TissueSpecific", method,norm.method, "seurat.PCA.dimension.selection.pdf" ))
  print(JackStrawPlot(tissue.raw.counts.object, dims = 1:15))
  print(ElbowPlot(tissue.raw.counts.object))
  dev.off()

  defined.pc = 1:15
  ########################################
  ## Cluster the cells (using pre-defined PC number)
  ########################################
  tissue.raw.counts.object <- FindNeighbors(tissue.raw.counts.object, reduction = "pca", dims = defined.pc, k.param=15, force.recalc = TRUE) #using pre-defined PCA number
  tissue.raw.counts.object <- FindClusters(tissue.raw.counts.object, resolution = 3, algorithm = 1) 
  head(Idents(tissue.raw.counts.object), 5)
  tissue.raw.counts.object$oriIdent = gsub("-.*","",colnames(tissue.raw.counts.object))

  ########################################
  ## Run non-linear dimensional reduction (UMAP/tSNE)  (using pre-defined PC number)
  ########################################
  Idents(tissue.raw.counts.object) = tissue.raw.counts.object$oriIdent
  ## Run UMAP
  tissue.raw.counts.object <- RunUMAP(tissue.raw.counts.object, dims = defined.pc)
  pdf(file.path(prewd, "Output", "TissueSpecific", method,norm.method, "seurat.cluster.umap.ori.ident.pdf" ), width = 10)
  print(DimPlot(tissue.raw.counts.object, reduction = "umap"))
  dev.off()

  ## Run Non-linear dimensional reduction (tSNE)
  tissue.raw.counts.object <- RunTSNE(object = tissue.raw.counts.object, dims = defined.pc)
  pdf(file.path(prewd, "Output", "TissueSpecific", method,norm.method, "seurat.cluster.tsne.ori.ident.pdf" ), width = 10)
  print(DimPlot(object = tissue.raw.counts.object, reduction = "tsne")) # note that you can set do.label=T to help label individual clusters
  dev.off()
  
  saveRDS(tissue.raw.counts.object, file = file.path(prewd, "ProcessedData","TissueObjects",method,norm.method, "tissue.raw.counts.object.seurat.rds"))
  
  # for(stage in gsub("_","-",stage.order)){
  build.stage.de.genes = function(stage){
    message(stage)
    #stage = "p25year"
    # dir.create(file.path(prewd, "Output", "TissueSpecific",method,norm.method,stage))
    tissue.raw.counts.object.selected= tissue.raw.counts.object[,grep(stage,colnames(tissue.raw.counts.object))]
    ##################################################################
    ## Finding differentially expressed genes (cluster biomarkers)
    ##################################################################
    # tissue.raw.counts.object.markers <- FindAllMarkers(tissue.raw.counts.object, only.pos = TRUE, test.use = "negbinom", logfc.threshold = 0.25) #DESeq2
    # tissue.raw.counts.object.markers.deseq2 <- FindAllMarkers(tissue.raw.counts.object.selected, only.pos = TRUE, test.use = "DESeq2", logfc.threshold = 0.25, min.cells.group=1) #
    tissue.raw.counts.object.markers.nb <- FindAllMarkers(tissue.raw.counts.object.selected, only.pos = TRUE, test.use = "negbinom", logfc.threshold = 0.25, min.cells.group=1, min.cells.feature=1) #
    tissue.raw.counts.object.markers = tissue.raw.counts.object.markers.nb #rbind(tissue.raw.counts.object.markers.deseq2, tissue.raw.counts.object.markers.nb[!rownames(tissue.raw.counts.object.markers.nb) %in% rownames(tissue.raw.counts.object.markers.deseq2),])
    
    promoter.tissue.specific.seurat.candidates = tissue.raw.counts.object.markers
    promoter.tissue.specific.seurat.object = tissue.raw.counts.object
    save(promoter.tissue.specific.seurat.candidates, file = file.path(prewd, "ProcessedData", "TissueObjects",method,norm.method, paste0("promoter.tissue.specific.seurat.candidates.", stage, ".RData") ))
    save(promoter.tissue.specific.seurat.object, file = file.path(prewd, "ProcessedData", "TissueObjects",method,norm.method, paste0("promoter.tissue.specific.seurat.object.", stage, ".RData") ))
    
    top10 <- tissue.raw.counts.object.markers %>% group_by(cluster) %>% top_n(n = 80, wt = avg_logFC)
    pdf(file.path(prewd, "Output", "TissueSpecific",method,norm.method, paste0("seurat.cluster.marker.doheatmap.ori.ident.", stage, ".pdf" )))
    print(DoHeatmap(tissue.raw.counts.object.selected, features = top10$gene) + NoLegend())
    dev.off()
  }
  
  registerDoParallel(cores = 5)
  foreach (stage = gsub("_","-",stage.order)) %dopar% build.stage.de.genes(stage)
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}


################################################
####Tissue specific analysis
################################################
###tissue specific identification
if(TRUE){
  # sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects",'human.sample.sheet.txt'),header=TRUE, stringsAsFactors = FALSE)
  all.stage = stage.order #unique(sample.info$stage) #table(sapply((strsplit(colnames(absolute.promoter.activity.star.mean),"_")),function(x){return(x[2])}))
  
  load(file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.RData'))
  
  ##label genes with only one promoter
  selected.promoters = promoterAnnotationData@promoterCoordinates[promoterAnnotationData@promoterCoordinates$promoterId %in% rownames(absolute.promoter.activity.star.mean)]
  selected.promoters.name = selected.promoters$geneId[duplicated(selected.promoters$geneId)]
  selected.promoters.twoplus = selected.promoters[selected.promoters$geneId %in% selected.promoters.name]
  
  ##In each tissue, compare with other tissue,  mean absolute promoter activity tissue/other>= 2; mean gene expression tissue/other<= 1.5
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'gene.expression.star.mean.RData'))
  gene.expression.star.mean = gene.expression.star.mean[,colnames(absolute.promoter.activity.star.mean)]
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'relative.promoter.activity.star.mean.RData'))
  relative.promoter.activity.star.mean[is.na(relative.promoter.activity.star.mean)] = 0
  relative.promoter.activity.star.mean = relative.promoter.activity.star.mean[,colnames(absolute.promoter.activity.star.mean)]
  
  for(stage in all.stage){
    #stage = "p25year" #"s9wpb_Male"
    load(file.path(prewd, "ProcessedData", "TissueObjects",method,norm.method, paste0("promoter.tissue.specific.seurat.candidates.", gsub("_","-",stage), ".RData") ))
    absolute.promoter.activity.star.mean.selected = absolute.promoter.activity.star.mean[,grep(stage,colnames(absolute.promoter.activity.star.mean))]
    relative.promoter.activity.star.mean.selected = relative.promoter.activity.star.mean[,grep(stage,colnames(relative.promoter.activity.star.mean))]
    gene.expression.star.mean.selected = gene.expression.star.mean[,grep(stage,colnames(gene.expression.star.mean))]
    gene.expression.star.mean.selected = gene.expression.star.mean.selected[!apply(gene.expression.star.mean.selected, 1, function(x)any(x==0)),]
    
    all.tissue = unique(sapply((strsplit(colnames(absolute.promoter.activity.star.mean.selected),"_")),function(x){return(x[1])}))
    promoter.tissue.specific.list=list()
    promoter.tissue.specific.granges=list()
    for(tissue in all.tissue){
      promoter.activity.star.mean.all = data.frame(absolute.tissue = rowMeans(absolute.promoter.activity.star.mean.selected[,grep(tissue,colnames(absolute.promoter.activity.star.mean.selected)),drop=FALSE]),
                                              absolute.other = rowMeans(absolute.promoter.activity.star.mean.selected[,-grep(tissue,colnames(absolute.promoter.activity.star.mean.selected)),drop=FALSE]),
                                              relative.tissue = rowMeans(relative.promoter.activity.star.mean.selected[,grep(tissue,colnames(relative.promoter.activity.star.mean.selected)),drop=FALSE]), 
                                              relative.other = rowMeans(relative.promoter.activity.star.mean.selected[,-grep(tissue,colnames(relative.promoter.activity.star.mean.selected)),drop=FALSE]), #apply(, 1, mean)
                                              # tau = absolute.promoter.activity.star.mean.tissue.specific.tau[rownames(absolute.promoter.activity.star.mean.selected),stage],
                                              # spm = absolute.promoter.activity.star.mean.tissue.specific.spm[rownames(absolute.promoter.activity.star.mean.selected),stage],
                                              stringsAsFactors = FALSE)
      
      gene.expression.star.mean.all = data.frame(expression.tissue = rowMeans(gene.expression.star.mean.selected[,grep(tissue,colnames(gene.expression.star.mean.selected)),drop=FALSE]),
                                            expression.other = rowMeans(gene.expression.star.mean.selected[,-grep(tissue,colnames(gene.expression.star.mean.selected)),drop=FALSE]),
                                            stringsAsFactors = FALSE)
      
      promoter.activity.star.mean.all.merge = merge(promoter.activity.star.mean.all,unique(promoterAnnotationData@promoterIdMapping[,c("promoterId","geneId")]),by.x=0,by.y="promoterId")
      promoter.activity.star.mean.all.merge = merge(promoter.activity.star.mean.all.merge,gene.expression.star.mean.all,by.x="geneId",by.y=0)
      rownames(promoter.activity.star.mean.all.merge) = promoter.activity.star.mean.all.merge$Row.names
      promoter.activity.star.mean.all.merge$absolute.divide = promoter.activity.star.mean.all.merge$absolute.tissue/promoter.activity.star.mean.all.merge$absolute.other
      promoter.activity.star.mean.all.merge$relative.subtract = promoter.activity.star.mean.all.merge$relative.tissue - promoter.activity.star.mean.all.merge$relative.other
      promoter.activity.star.mean.all.merge$expression.divide = promoter.activity.star.mean.all.merge$expression.tissue/promoter.activity.star.mean.all.merge$expression.other
      promoter.activity.star.mean.all.merge[is.na(promoter.activity.star.mean.all.merge)] = 0
      promoter.activity.star.mean.all.merge = promoter.activity.star.mean.all.merge[rownames(promoter.activity.star.mean.all.merge) %in% selected.promoters.twoplus$promoterId,]
      
      promoter.activity.tissue.specific = merge(promoter.activity.star.mean.all.merge[,-2], promoter.tissue.specific.seurat.candidates[promoter.tissue.specific.seurat.candidates$cluster == tissue, ], by = 0)
      rownames(promoter.activity.tissue.specific) = promoter.activity.tissue.specific$Row.names
      
      promoter.activity.tissue.specific.intersect = promoter.activity.tissue.specific[promoter.activity.tissue.specific$absolute.tissue>=0.25 & 
                                                                             # promoter.activity.tissue.specific$absolute.other>=0.25 &
                                                                             promoter.activity.tissue.specific$relative.tissue>=0.2 &
                                                                             # promoter.activity.tissue.specific$relative.other>=0.25 & 
                                                                             promoter.activity.tissue.specific$expression.tissue>=0.75 &
                                                                             promoter.activity.tissue.specific$expression.other>=0.75 &
                                                                             # promoter.activity.tissue.specific$absolute.divide>=2 &
                                                                            promoter.activity.tissue.specific$expression.divide<=3 &
                                                                             promoter.activity.tissue.specific$relative.subtract>=0.2,
                                                                            ,drop=FALSE] #(promoter.activity.tissue.specific$tau >=0.3 | promoter.activity.tissue.specific$spm >=0.3)
      # 
      promoterAnnotationDataGranges = promoterAnnotationData@promoterCoordinates
      promoterAnnotationDataGranges.selected = promoterAnnotationDataGranges[promoterAnnotationDataGranges$promoterId %in% rownames(promoter.activity.tissue.specific.intersect)]
      # absolute.promoter.activity.star.mean.selected[paste0("prmtr.2516",c(1,3,5,6,7)),]
      # dcounts.selected[paste0("prmtr.2516",c(1,3,5,6,7)),]
      
      promoter.tissue.specific.list[[tissue]] = promoter.activity.tissue.specific.intersect
      promoter.tissue.specific.granges[[tissue]] = promoterAnnotationDataGranges.selected
    }
    save(promoter.tissue.specific.list,file=file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste0("promoter.tissue.specific.list.", stage, ".RData") ))
    save(promoter.tissue.specific.granges,file=file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste0("promoter.tissue.specific.granges.", stage, ".RData")) )
    
    
  }
  # save(promoter.tissue.specific.list,file=file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,'promoter.tissue.specific.list.RData'))
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

##TS example
if(TRUE){
  for(stage in stage.order){
    # stage="p25year"
    load(file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
    load(file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste0("promoter.tissue.specific.granges.", stage, ".RData")) )
    load(file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste0("promoter.tissue.specific.list.", stage, ".RData")) )
    for (name in names(promoter.tissue.specific.granges)){
      tmp = promoter.tissue.specific.granges[[name]]
      tmp$name = tmp$promoterId
      rtracklayer::export(tmp, con=file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste(stage, name, "bed", sep=".")) )
    }
  }
  # promoter.anno = promoterAnnotationData@promoterCoordinates
  # names(promoter.anno) = promoter.anno$promoterId
  # brain.specific = promoter.tissue.specific.list[[grep("Brain",names(promoter.tissue.specific.list))]]
  # export(promoter.anno[rownames(brain.specific[brain.specific$expression.divide<=1.5,])] + 1000,con="Mammalian2019Promoter/CM2019RNAseqHuman/PygenomeSnapshot/brain.specific.bed")
  rm(list=setdiff(ls(),basic.objects))
  gc()
}


####################################################################################
###statistics of tissue specific gene number in each developmental stage
####################################################################################
if(TRUE){
  promoter.tissue.specific.number.statistics = rbind()
  for (stage in stage.order){
    load(file=file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste0("promoter.tissue.specific.list.", stage, ".RData") ))
    promoter.tissue.specific.number=sapply(promoter.tissue.specific.list, dim)
    promoter.tissue.specific.number.per.stage=data.frame(tissue = colnames(promoter.tissue.specific.number), value = promoter.tissue.specific.number[1, ], stage=stage)
    rownames(promoter.tissue.specific.number.per.stage) = NULL
    promoter.tissue.specific.number.statistics = rbind(promoter.tissue.specific.number.statistics, promoter.tissue.specific.number.per.stage)
  }

  ##barplot show number distribution in each tissue
  p1<-ggplot(promoter.tissue.specific.number.statistics, aes(x=stage,y=value,fill=stage)) + 
    geom_bar(position="dodge",stat="identity") + 
    # scale_x_discrete(breaks = stage.order) +
    labs(y = "Number of alternative promoters") + 
    facet_grid(tissue ~ .) +
    theme_bw() + theme(axis.text.x=element_text(angle = 90),
                       panel.grid.major.y = element_blank(),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       legend.position = "none")
  p1
  ggsave(p1,file=file.path(prewd,"Output","TissueSpecific",method,norm.method,'promoter.tissue.specific.number.statistics.barplot.pdf'))
  
  
  pline=ggplot(promoter.tissue.specific.number.statistics, aes(x=stage,y=value, group=tissue, shape=tissue,color=tissue)) + 
    geom_line(stat = "smooth", method = "loess",size=1.5) + 
    geom_point(size=1, fill="white") +
    scale_shape_manual(values=rep(19,length(unique(promoter.tissue.specific.number.statistics$tissue))))+
    # scale_x_discrete(labels=stage.order) +
    scale_fill_manual(values =c("#d59bf6", "#c264fe", "#a82ffc", "#b61aae", "#590d82", "#6a2c70", "#521262")) + 
    labs(y = "Number of alternative promoters")+ 
    theme_classic()+theme(axis.text.x=element_text(angle = 90))+
    ylim(-12,155)
  
  pline
  ggsave(pline,file=file.path(prewd,"Output","TissueSpecific",method,norm.method,'promoter.tissue.specific.number.statistics.lineplot.pdf'))
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}


################################################
###heatmap of tissue specific promoters across stages and tissues
################################################
if(TRUE){
  promoter.tissue.specific.list.all = list()
  for (stage in stage.order){
    load(file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste0("promoter.tissue.specific.list.", stage, ".RData") ))
    names(promoter.tissue.specific.list) = paste(names(promoter.tissue.specific.list),stage,sep="_")
    promoter.tissue.specific.list.all = c(promoter.tissue.specific.list.all, promoter.tissue.specific.list)
  }
  promoter.tissue.specific.list = promoter.tissue.specific.list.all
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.RData'))
  # sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects",'human.sample.sheet.txt'),header=TRUE, stringsAsFactors = FALSE)
  
  #heatmap for all, arragned based on tissue
  promoter.tissue.specific.allpromoter= melt(sapply((promoter.tissue.specific.list),function(x){return(rownames(x))}))
  promoter.tissue.specific.allpromoter.activity= absolute.promoter.activity.star.mean[as.character(promoter.tissue.specific.allpromoter$value), sort(colnames(absolute.promoter.activity.star.mean))]
  pheatmap(promoter.tissue.specific.allpromoter.activity, 
           color = colorRampPalette(c("#a1e6e3","white","#f64b3c"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE,show_colnames = T,scale = "row", width = 20, height=20, legend = TRUE,
           filename = file.path(prewd,"Output","TissueSpecific",method,norm.method,"promoter.tissue.specific.heatmap.all.tissue.pdf"))
  
  #heatmap for all, arragned based on stage
  promoter.tissue.specific.allpromoter.activity.stage = cbind()
  for(sample in names(promoter.tissue.specific.list)){
    promoter.tissue.specific.stage.allpromoter = promoter.tissue.specific.allpromoter.activity[ ,sample] 
    promoter.tissue.specific.allpromoter.activity.stage = cbind(promoter.tissue.specific.allpromoter.activity.stage, promoter.tissue.specific.stage.allpromoter)
  }
  colnames(promoter.tissue.specific.allpromoter.activity.stage) = names(promoter.tissue.specific.list)
  
  pheatmap(promoter.tissue.specific.allpromoter.activity.stage, 
           color = colorRampPalette(c("#a1e6e3","white","#f64b3c"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE,show_colnames = T,scale = "row", width = 20, height=20, legend = TRUE,
           filename = file.path(prewd,"Output","TissueSpecific",method,norm.method,"promoter.tissue.specific.heatmap.all.stage.pdf"))
  
  #heatmap for each stage (b4week)
  for (stage in stage.order){
    # stage="b4week"
    promoter.tissue.specific.stage=promoter.tissue.specific.list[grep(stage,names(promoter.tissue.specific.list))]
    promoter.tissue.specific.stage.allpromoter=melt(sapply((promoter.tissue.specific.stage),function(x){return(rownames(x))}))
    promoter.tissue.specific.stage.allpromoter.activity=absolute.promoter.activity.star.mean[,grep(stage,colnames(absolute.promoter.activity.star.mean))]
    promoter.tissue.specific.stage.allpromoter.activity=promoter.tissue.specific.stage.allpromoter.activity[as.character(promoter.tissue.specific.stage.allpromoter$value), ]
    
    pheatmap(promoter.tissue.specific.stage.allpromoter.activity, 
             color = colorRampPalette(c("#a1e6e3","white","#f64b3c"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
             cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=8, legend = TRUE,
             filename = file.path(prewd,"Output","TissueSpecific",method,norm.method,paste0("promoter.tissue.specific.heatmap.", stage, ".pdf")))
  }
  
  #heatmap for each tissue (Brain)
  for (tissue in tissue.order){
    # tissue="Brain"
    promoter.tissue.specific.tissue=promoter.tissue.specific.list[grep(tissue,names(promoter.tissue.specific.list))]
    promoter.tissue.specific.tissue.allpromoter=melt(sapply((promoter.tissue.specific.tissue),function(x){return(rownames(x))}))
    promoter.tissue.specific.tissue.allpromoter.activity=absolute.promoter.activity.star.mean[,grep(tissue,colnames(absolute.promoter.activity.star.mean))]
    promoter.tissue.specific.tissue.allpromoter.activity=promoter.tissue.specific.tissue.allpromoter.activity[as.character(promoter.tissue.specific.tissue.allpromoter$value), ]
    
    pheatmap(promoter.tissue.specific.tissue.allpromoter.activity, 
             color = colorRampPalette(c("#a1e6e3","white","#f64b3c"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
             cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=8, legend = TRUE,
             filename = file.path(prewd,"Output","TissueSpecific",method,norm.method,paste0("promoter.tissue.specific.heatmap.", tissue, ".pdf")))
  }
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}



###########################################################
###heatmap of shared gene number between different stages in one tissue
###########################################################
if(TRUE){
  # load AbsolutePromoterActivity and sample info
  # sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects",'human.sample.sheet.txt'),header=TRUE, stringsAsFactors = FALSE)
  promoter.tissue.specific.list.all = list()
  for (stage in stage.order){
    load(file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste0("promoter.tissue.specific.list.", stage, ".RData") ))
    names(promoter.tissue.specific.list) = paste(names(promoter.tissue.specific.list),stage,sep="_")
    promoter.tissue.specific.list.all = c(promoter.tissue.specific.list.all, promoter.tissue.specific.list)
  }
  promoter.tissue.specific.list = promoter.tissue.specific.list.all
  
  for (tissue in tissue.order){
    # tissue = "Brain"
    tissue.specific=promoter.tissue.specific.list[grep(tissue,names(promoter.tissue.specific.list))]
    tissue.overlap.allpromoter=sapply((tissue.specific),function(x){return(rownames(x))})
    tissue.overlap.allpromoter.heatmap=sapply(seq_len(length(tissue.overlap.allpromoter)), function(x) 
      sapply(seq_len(length(tissue.overlap.allpromoter)), function(y) length(intersect(unlist(tissue.overlap.allpromoter[x]), unlist(tissue.overlap.allpromoter[y])))))
    # tissue.name=unique(sample.info[grep(tissue, sample.info$condition), "condition"])
    colnames(tissue.overlap.allpromoter.heatmap)=names(tissue.specific)
    rownames(tissue.overlap.allpromoter.heatmap)=names(tissue.specific)
    tissue.overlap.allpromoter.heatmap=log2(tissue.overlap.allpromoter.heatmap)
    # Melt the correlation matrix
    melted_cormat <- melt(tissue.overlap.allpromoter.heatmap, na.rm = TRUE)
    melted_cormat[melted_cormat == -Inf] = 0
    # Heatmap
    p1=ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+ 
      geom_vline(xintercept = 14.5, linetype="dotted")+
      geom_hline(yintercept = 14.5, linetype="dotted")+
      scale_fill_gradient2(low = "white", high = "red",  
                           space = "Lab", 
                           name="Pearson\nCorrelation") +
      theme_minimal()+ theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.3,size = 5, hjust = 0.3),
            axis.text.y = element_text(size = 5, hjust = 0.3)) +
      labs(x="", y = "") +
      coord_fixed() 
    p1
    ggsave(p1,file=file.path(prewd,"Output","TissueSpecific",method,norm.method,paste0("overlap.of.tissue.specific.promoters.across.stages.", tissue, ".pdf")))
  }
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}


################################################
###example Testis_p15year major/minor/expression (Aamp, Tuba4a, mir8114,Rfxap)
################################################
if(TRUE){
  # load RelativePromoterActivity
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'relativePromoterActivity.star.RData'))
  # load GeneExpression
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'geneExpression.star.RData'))
  # load sample info
  # sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects",'human.sample.sheet.txt'),header=TRUE, stringsAsFactors = FALSE)
  # load AbsolutePromoterActivity and sample info
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolutePromoterActivity.star.RData'))
  promoter.tissue.specific.list.all = list()
  for (stage in stage.order){
    load(file.path(prewd,"ProcessedData","TissueObjects",method,norm.method,paste0("promoter.tissue.specific.list.", stage, ".RData") ))
    names(promoter.tissue.specific.list) = paste(names(promoter.tissue.specific.list),stage,sep="_")
    promoter.tissue.specific.list.all = c(promoter.tissue.specific.list.all, promoter.tissue.specific.list)
  }
  promoter.tissue.specific.list = promoter.tissue.specific.list.all
  
  stage = "p15year" #Testis_ #Liver_b4week
  tissue = "Testis"
  major.pro = "prmtr.14617" # "prmtr.37697"
  minor.pro = "prmtr.14616" # "prmtr.37694"
  gene.name = "ENSG00000114573"
  
  absolutePromoterActivity.star.major=melt(absolutePromoterActivity.star[major.pro,grep(stage,colnames(absolutePromoterActivity.star))])
  absolutePromoterActivity.star.major$variable = gsub(paste0(tissue,"_.*"),tissue,absolutePromoterActivity.star.major$variable)
  absolutePromoterActivity.star.major$variable = gsub(".*_.*","OtherTissues",absolutePromoterActivity.star.major$variable)
  
  p1<-ggplot(absolutePromoterActivity.star.major, aes(x=variable, y=value, color=variable,fill=variable)) + 
    geom_boxplot(alpha=1,outlier.shape = NA,width=0.4) +
    scale_fill_manual(values =c("#ff1e56", "#ffac41")) + 
    scale_color_manual(values =rep("black",2))+
    labs(x="", y = "Promoter Activity")+
    stat_compare_means(label = "p.signif",method = "wilcox.test",ref.group = tissue,na.rm = TRUE,size=rel(8)) + 
    ggtitle(major.pro)+theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
  p1
  
  ggsave(p1,file=file.path(prewd,"Output","TissueSpecific",method,norm.method,paste("promoter.activity", tissue, stage, major.pro, "pdf", sep=".")),width=4,height = 4,units = "in",dpi=300)
  
  #look for minor promoter prmtr.37694
  absolutePromoterActivity.star.minor=melt(absolutePromoterActivity.star[minor.pro,grep(stage,colnames(absolutePromoterActivity.star))])
  absolutePromoterActivity.star.minor$variable = gsub(paste0(tissue,"_.*"),tissue,absolutePromoterActivity.star.minor$variable)
  absolutePromoterActivity.star.minor$variable = gsub(".*_.*","OtherTissues",absolutePromoterActivity.star.minor$variable)
  
  p2<-ggplot(absolutePromoterActivity.star.minor, aes(x=variable, y=value, color=variable,fill=variable)) + 
    geom_boxplot(alpha=1,outlier.shape = NA,width=0.4) +
    scale_fill_manual(values =c("#ff1e56", "#ffac41")) + 
    scale_color_manual(values =rep("black",2))+
    labs(x="", y = "Promoter Activity")+
    stat_compare_means(label = "p.signif",method = "wilcox.test",ref.group = tissue,na.rm = TRUE,size=rel(8)) + 
    ggtitle(minor.pro)+theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
  p2
  
  ggsave(p2,file=file.path(prewd,"Output","TissueSpecific",method,norm.method,paste("promoter.activity", tissue, stage, minor.pro, "pdf", sep=".")),width=4,height = 4,units = "in",dpi=300)
  
  
  #gene expression
  
  geneExpression.star.selected=melt(geneExpression.star[gene.name,grep(stage,colnames(geneExpression.star))])
  geneExpression.star.selected$variable = gsub(paste0(tissue,"_.*"),tissue,geneExpression.star.selected$variable)
  geneExpression.star.selected$variable = gsub(".*_.*","OtherTissues",geneExpression.star.selected$variable)
  
  p3<-ggplot(geneExpression.star.selected, aes(x=variable, y=value, color=variable,fill=variable)) + 
    geom_boxplot(alpha=1,outlier.shape = NA,width=0.4) +
    scale_fill_manual(values =c("#ff1e56", "#ffac41")) + 
    scale_color_manual(values =rep("black",2))+
    labs(x="", y = "Gene expression")+
    stat_compare_means(label = "p.signif",method = "wilcox.test",ref.group = tissue,na.rm = TRUE,size=rel(8)) + 
    ggtitle(gene.name)+theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
  p3
  
  ggsave(p3,file=file.path(prewd,"Output","TissueSpecific",method,norm.method,paste("promoter.activity", tissue, stage, gene.name, "pdf", sep=".")),width=4,height = 4,units = "in",dpi=300)
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}



# #########tissue specific plot
# load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'promoter.tissue.specific.list.RData'))
# promoter.tissue.specific=sapply(promoter.tissue.specific.list, dim)
# promoter.tissue.specific=t(promoter.tissue.specific[1, ])
# promoter.tissue.specific=melt(promoter.tissue.specific)
# promoter.tissue.specific=promoter.tissue.specific[, -1]
# brain.plot = promoter.tissue.specific[grep("Brain",promoter.tissue.specific$Var2), ]
# Cerebellum.plot = promoter.tissue.specific[grep("Cerebellum",promoter.tissue.specific$Var2), ]
# Heart.plot = promoter.tissue.specific[grep("Heart",promoter.tissue.specific$Var2), ]
# Kidney.plot = promoter.tissue.specific[grep("Kidney",promoter.tissue.specific$Var2), ]
# Ovary.plot = promoter.tissue.specific[grep("Ovary",promoter.tissue.specific$Var2), ]
# Testis.plot = promoter.tissue.specific[grep("Testis",promoter.tissue.specific$Var2), ]
# Liver.plot = promoter.tissue.specific[grep("Liver",promoter.tissue.specific$Var2), ]


##Tissue specific score, Not used anymore
##Tau method
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
fTau <- function(x)
{
  if(all(!is.na(x)))
  {if(min(x, na.rm=TRUE) >= 0)
  {if(max(x)!=0)
  {x <- (1-(x/max(x)))
  res <- sum(x, na.rm=TRUE)
  res <- res/(length(x)-1)
  } else {res <- 0}
  } else {res <- NA } #print("Expression values have to be positive!")
  } else {res <- NA } #print("No data for this gene avalable.")
  return(res)
}
###***###***###
##SPM method: 
###+++###
#SPM score from TISGED
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fSpm <- function(x)
{if(all(!is.na(x)))
{if(min(x, na.rm=TRUE) >= 0)
{if(sum(x) !=0)
{spm <- x^2/(x%*%x)
res <- max(spm) #Modification:To bring to normalized scale. Choose max
} else {res <- 0}	 		
} else {res <- NA} #print("Expression values have to be positive!")
} else {res <- NA} #print("No data for this gene avalable.")
  return(res)
}
# ###tissue specific identification
# if(TRUE){
#   sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects",'human.sample.sheet.txt'),header=TRUE, stringsAsFactors = FALSE)
#   load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.RData'))
#   
#   ###tissue specific score calculation in each developmental stage using tau and spm###
#   all.stage = unique(sample.info$stage) #table(sapply((strsplit(colnames(absolute.promoter.activity.star.mean),"_")),function(x){return(x[2])}))
#   absolute.promoter.activity.star.mean.tissue.specific.tau = cbind()
#   for(stage in all.stage){
#     absolute.promoter.activity.star.mean.tissue.specific.tau = cbind(absolute.promoter.activity.star.mean.tissue.specific.tau, apply(absolute.promoter.activity.star.mean[,grep(stage,colnames(absolute.promoter.activity.star.mean))], 1, fTau))
#   }
#   colnames(absolute.promoter.activity.star.mean.tissue.specific.tau) = all.stage
#   save(absolute.promoter.activity.star.mean.tissue.specific.tau,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.tissue.specific.tau.RData'))
#   
#   absolute.promoter.activity.star.mean.tissue.specific.spm = cbind()
#   for(stage in all.stage){
#     absolute.promoter.activity.star.mean.tissue.specific.spm = cbind(absolute.promoter.activity.star.mean.tissue.specific.spm, apply(absolute.promoter.activity.star.mean[,grep(stage,colnames(absolute.promoter.activity.star.mean))], 1, fSpm))
#   }
#   colnames(absolute.promoter.activity.star.mean.tissue.specific.spm) = all.stage
#   save(absolute.promoter.activity.star.mean.tissue.specific.spm,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.tissue.specific.spm.RData'))
#   
#   ##In each tissue, compare with other tissue,  mean absolute promoter activity tissue/other>= 2; mean gene expression tissue/other<= 1.5
#   load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'gene.expression.star.mean.RData'))
#   gene.expression.star.mean = gene.expression.star.mean[,colnames(absolute.promoter.activity.star.mean)]
#   load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'relative.promoter.activity.star.mean.RData'))
#   relative.promoter.activity.star.mean[is.na(relative.promoter.activity.star.mean)] = 0
#   relative.promoter.activity.star.mean = relative.promoter.activity.star.mean[,colnames(absolute.promoter.activity.star.mean)]
#   
#   load(file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
#   load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.tissue.specific.tau.RData'))
#   load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'absolute.promoter.activity.star.mean.tissue.specific.spm.RData'))
#   
#   promoter.tissue.specific.list=list()
#   for(stage in all.stage){
#     absolute.promoter.activity.star.mean.selected = absolute.promoter.activity.star.mean[,grep(stage,colnames(absolute.promoter.activity.star.mean))]
#     relative.promoter.activity.star.mean.selected = relative.promoter.activity.star.mean[,grep(stage,colnames(relative.promoter.activity.star.mean))]
#     gene.expression.star.mean.selected = gene.expression.star.mean[,grep(stage,colnames(gene.expression.star.mean))]
#     
#     all.tissue = unique(sapply((strsplit(colnames(absolute.promoter.activity.star.mean.selected),"_")),function(x){return(x[1])}))
#     for(tissue in all.tissue){
#       promoter.activity.star.mean.all = data.frame(absolute.tissue = rowMeans(absolute.promoter.activity.star.mean.selected[,grep(tissue,colnames(absolute.promoter.activity.star.mean.selected)),drop=FALSE]),
#                                               absolute.other = rowMeans(absolute.promoter.activity.star.mean.selected[,-grep(tissue,colnames(absolute.promoter.activity.star.mean.selected)),drop=FALSE]),
#                                               relative.tissue = rowMeans(relative.promoter.activity.star.mean.selected[,grep(tissue,colnames(relative.promoter.activity.star.mean.selected)),drop=FALSE]),
#                                               relative.other = rowMeans(relative.promoter.activity.star.mean.selected[,-grep(tissue,colnames(relative.promoter.activity.star.mean.selected)),drop=FALSE]),
#                                               tau = absolute.promoter.activity.star.mean.tissue.specific.tau[rownames(absolute.promoter.activity.star.mean.selected),stage],
#                                               spm = absolute.promoter.activity.star.mean.tissue.specific.spm[rownames(absolute.promoter.activity.star.mean.selected),stage],
#                                               stringsAsFactors = FALSE)
#       
#       gene.expression.star.mean.all = data.frame(expression.tissue = rowMeans(gene.expression.star.mean.selected[,grep(tissue,colnames(gene.expression.star.mean.selected)),drop=FALSE]),
#                                             expression.other = rowMeans(gene.expression.star.mean.selected[,-grep(tissue,colnames(gene.expression.star.mean.selected)),drop=FALSE]),
#                                             stringsAsFactors = FALSE)
#       
#       promoter.activity.star.mean.all.merge = merge(promoter.activity.star.mean.all,unique(promoterAnnotationData@promoterIdMapping[,c("promoterId","geneId")]),by.x=0,by.y="promoterId")
#       promoter.activity.star.mean.all.merge = merge(promoter.activity.star.mean.all.merge,gene.expression.star.mean.all,by.x="geneId",by.y=0)
#       rownames(promoter.activity.star.mean.all.merge) = promoter.activity.star.mean.all.merge$Row.names
#       promoter.activity.star.mean.all.merge$absolute.divide = promoter.activity.star.mean.all.merge$absolute.tissue/promoter.activity.star.mean.all.merge$absolute.other
#       promoter.activity.star.mean.all.merge$expression.divide = promoter.activity.star.mean.all.merge$expression.tissue/promoter.activity.star.mean.all.merge$expression.other
#       promoter.activity.star.mean.all.merge[is.na(promoter.activity.star.mean.all.merge)] = 0
#       
#       promoter.activity.tissue.specific = promoter.activity.star.mean.all.merge[promoter.activity.star.mean.all.merge$absolute.tissue>=0.25 & 
#                                                                              promoter.activity.star.mean.all.merge$absolute.other>=0.25 &
#                                                                              promoter.activity.star.mean.all.merge$relative.tissue>=0.25 & 
#                                                                              promoter.activity.star.mean.all.merge$relative.other>=0.25 & 
#                                                                              promoter.activity.star.mean.all.merge$expression.tissue>=0.75 & 
#                                                                              promoter.activity.star.mean.all.merge$expression.other>=0.75 & 
#                                                                              promoter.activity.star.mean.all.merge$absolute.divide>=2 & 
#                                                                              promoter.activity.star.mean.all.merge$expression.divide<2 & 
#                                                                              (promoter.activity.star.mean.all.merge$tau >=0.3 | promoter.activity.star.mean.all.merge$spm >=0.3),,drop=FALSE]
#       
#       promoter.tissue.specific.list[[paste(tissue,stage,sep="_")]] = promoter.activity.tissue.specific
#     }
#     
#   }
#   save(promoter.tissue.specific.list,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,'promoter.tissue.specific.list.RData'))
#   
#   
#   rm(list=setdiff(ls(),basic.objects))
#   gc()
# }
# 
################################################
###Tissue specifc analysis with Seurat per stage
################################################

# if(TRUE){
#   # ##using DEXseq counts
#   # load(file.path(prewd, "ProcessedData", "PromoterObjects",method,norm.method, "dexseq.counts.dcounts.RData"))
#   # load(file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.RData"))
#   # load(file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.transcripts.RData"))
#   # load(file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.promoter.correspond.RData"))
#   # dcounts.change.name = merge(dexseq.exoninfo.promoter.correspond, dcounts, by.x = "dexseq", by.y=0)  #take only promoters
#   # rownames(dcounts.change.name) = dcounts.change.name$promoterId
#   # # dcounts.change.name = dcounts.change.name[dcounts.change.name$internalPromoter==FALSE,]
#   # dcounts = dcounts.change.name[,-c(1:4)]
#   
#   ##load counts
#   load(file.path(prewd,"ProcessedData","PromoterObjects",method,'promoterCounts.star.RData'))
#   dcounts = promoterCounts.star
#   
#   dir.create(file.path(prewd, "ProcessedData", "TissueObjects"))
#   dir.create(file.path(prewd, "Output", "TissueSpecific",method))
#   dir.create(file.path(prewd, "Output", "TissueSpecific",method,norm.method))
#   
#   
#   for (stage in stage.order){ #unique(unlist(sapply(strsplit(colnames(dcounts), "_"), function(x){return(x[2])})))
#     ##Each developmental stage run seperately
#     # stage = "p5year" #b4week
#     dir.create(file.path(prewd, "Output", "TissueSpecific",method,norm.method, stage))
#     
#     ##Extract stage sample and build sample matrix
#     dcounts.selected = dcounts[,grep(stage, colnames(dcounts))]
#     # dcounts.selected = dcounts.selected[rownames(dcounts.selected) %in% dexseq.exoninfo.promoter.correspond$dexseq,]  #take only promoters
#     # dcounts.selected = dcounts.selected[rowSums(dcounts.selected) != 0 , ]
#     
#     tissue.raw.counts.object = CreateSeuratObject(counts = dcounts.selected, project = "HumanTissue", min.cells = 1, min.features = 0)
#     
#     ##########################
#     ##Quality Control
#     ##########################
#     # Visualize QC metrics as a violin plot
#     pdf(file.path(prewd, "Output", "TissueSpecific",method,norm.method, stage, "seurat.feature.vlnplot.pdf" ))
#     print(VlnPlot(tissue.raw.counts.object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) )
#     dev.off()
#     
#     # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#     # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#     plot2 <- FeatureScatter(tissue.raw.counts.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#     pdf(file.path(prewd, "Output", "TissueSpecific",method,norm.method, stage, "seurat.feature.scatterplot.pdf" ), width=14, height = 7)
#     print(plot2)
#     dev.off()
#     
#     ##subselect the data
#     # tissue.raw.counts.object <- subset(tissue.raw.counts.object, subset = (nFeature_RNA > min.feature & nFeature_RNA < max.feature & percent.mt < percent.mt) ) #nFeature_RNA < 2500 & nCount_RNA < 300 &
#     # expr <- FetchData(object = tissue.raw.counts.object, vars = "nFeature_RNA")
#     # tissue.raw.counts.object = tissue.raw.counts.object[, which(x = expr > min.feature & expr < max.feature)]
#     
#     ##########################
#     ##Normalization
#     ##########################
#     tissue.raw.counts.object <- NormalizeData(tissue.raw.counts.object, normalization.method = "LogNormalize", scale.factor = 10000)
#     
#     ##########################
#     ##feature selection (highly variable features)
#     ##########################
#     tissue.raw.counts.object <- FindVariableFeatures(tissue.raw.counts.object, selection.method = "vst", nfeatures = 2000)
#     
#     # plot variable features with and without labels
#     plot1 <- VariableFeaturePlot(tissue.raw.counts.object)
#     plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(tissue.raw.counts.object), 20), repel = TRUE)
#     pdf(file.path(prewd, "Output", "TissueSpecific",method,norm.method, stage, "seurat.feature.selection.pdf" ))
#     print(plot2)
#     dev.off()
#     
#     ##########################################################################
#     ##Scaling the data and removing unwanted sources of variation
#     ##########################################################################
#     tissue.raw.counts.object <- ScaleData(object = tissue.raw.counts.object, features = rownames(tissue.raw.counts.object), vars.to.regress = c("nCount_RNA"))
#     
#     ##################################################################
#     ## Finding differentially expressed genes (cluster biomarkers)
#     ##################################################################
#     # tissue.raw.counts.object.markers <- FindAllMarkers(tissue.raw.counts.object, only.pos = TRUE, test.use = "negbinom", logfc.threshold = 0.25) #DESeq2
#     tissue.raw.counts.object.markers.deseq2 <- FindAllMarkers(tissue.raw.counts.object, only.pos = TRUE, test.use = "DESeq2", logfc.threshold = 1, min.cells.group=1) #
#     tissue.raw.counts.object.markers.poisson <- FindAllMarkers(tissue.raw.counts.object, only.pos = TRUE, test.use = "poisson", logfc.threshold = 1, min.cells.group=1) #
#     tissue.raw.counts.object.markers = rbind(tissue.raw.counts.object.markers.deseq2, tissue.raw.counts.object.markers.poisson[!rownames(tissue.raw.counts.object.markers.poisson) %in% rownames(tissue.raw.counts.object.markers.deseq2),])
#     tissue.raw.counts.object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) ##The first 2 for each cluster
#     
#     promoter.tissue.specific.seurat.candidates = tissue.raw.counts.object.markers
#     promoter.tissue.specific.seurat.object = tissue.raw.counts.object
#     save(promoter.tissue.specific.seurat.candidates, file = file.path(prewd, "ProcessedData", "TissueObjects",paste0("promoter.tissue.specific.seurat.candidates.", stage, ".RData") ))
#     save(promoter.tissue.specific.seurat.object, file = file.path(prewd, "ProcessedData", "TissueObjects",paste0("promoter.tissue.specific.seurat.object.", stage, ".RData") ))
#     
#     top10 <- tissue.raw.counts.object.markers %>% group_by(cluster) %>% top_n(n = 80, wt = avg_logFC)
#     pdf(file.path(prewd, "Output", "TissueSpecific",method,norm.method, stage, "seurat.cluster.marker.doheatmap.ori.ident.pdf" ))
#     print(DoHeatmap(tissue.raw.counts.object, features = top10$gene) + NoLegend())
#     dev.off()
#   }
#   
#   rm(list=setdiff(ls(),basic.objects))
#   gc()
# }


