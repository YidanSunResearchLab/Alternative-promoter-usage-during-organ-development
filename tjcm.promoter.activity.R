#Fig. 1b-f and Fig. S1a

# prewd="Mammalian2019Promoter/CM2019RNAseqHuman"
prewd="Mammalian2019Promoter/CM2019RNAseqMouse"
mc.cores= 20
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
library(ggpubr)

method="Proactiv"
# method="Dexseq"
# norm.method = "edger"
norm.method = "deseq2"

normalizePromoterReadCounts = function (promoterReadCounts, sizefactor=NULL) {
  if (ncol(promoterReadCounts) == 1) {return(promoterReadCounts)}
  activePromoters <- which(!is.na(promoterReadCounts[, 1]))
  colData <- data.frame(sampleLabels = colnames(promoterReadCounts))
  rownames(colData) <- colnames(promoterReadCounts)
  print("Calculating normalized read counts...")
  if (requireNamespace("DESeq2", quietly = TRUE) == FALSE) {
    stop("Error: DESeq2 is not installed! For normalization DESeq2 is needed, please install it.")
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = promoterReadCounts[activePromoters,], colData = colData, design = ~1)
  dds <- DESeq2::estimateSizeFactors(dds)
  if(!is.null(sizefactor)){DESeq2::sizeFactors(dds) = 1/sizefactor; message("sizefactor from outside")}
  promoterReadCounts[activePromoters, ] <- DESeq2::counts(dds,normalized = TRUE)
  return(promoterReadCounts)
}

theme_classic.add = theme_classic()+ theme(legend.position = "none",strip.background = element_blank(),strip.text = element_blank(), plot.title = element_text(hjust = 0.5, family = "Helvetica"),axis.text = element_text(colour = "black", family = "Helvetica"), axis.title=element_text(colour = "black", family = "Helvetica"), text = element_text(colour = "black", family = "Helvetica")) #axis.text = element_text(colour = "black"), axis.title=element_text(colour = "black"), 
basic.objects=c("basic.objects",ls())

################################################
#####Estimate Promoter Activity (STAR alignment) with DEseq2 size factor
################################################
if(TRUE){
  
  dir.create(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method))
  
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,"promoterCounts.star.RData"))
  load(file.path(prewd,"ProcessedData","AnnotationObjects","promoterAnnotationData.RData"))
  sizefactor = read.delim(file.path(prewd,"bamCoverage",paste0(norm.method,".dexseq.counts.txt")), row.names = 1, sep=":", stringsAsFactors = FALSE)
  rownames(sizefactor) = sample.info[rownames(sizefactor), "newnamegender"]
  
  if(norm.method == "deseq2"){
    # Normalize promoter read counts by DESeq2
    normalizedPromoterCounts.star <- normalizePromoterReadCounts(promoterCounts.star, 
                                                                 sizefactor=sizefactor[colnames(promoterCounts.star),1])
    save(normalizedPromoterCounts.star,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"normalizedPromoterCounts.star.RData"))
  } else if(norm.method == "edger"){
    # Normalize promoter read counts by edgeR
    library(edgeR)
    cds = DGEList(counts=as.matrix(promoterCounts.star[which(!is.na(promoterCounts.star[, 1])),]),group=factor(1:ncol(promoterCounts.star)))
    cds = calcNormFactors(cds, method=c("TMM"))
    cds$samples$norm.factors = 1/sizefactor[rownames(cds$samples),1]
    normalizedPromoterCounts.star = cpm(cds, normalized.lib.sizes = TRUE)
    save(normalizedPromoterCounts.star, file = file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"normalizedPromoterCounts.star.RData"))
  }
  
  
  # Calculate absolute promoter activity
  absolutePromoterActivity.star <- getAbsolutePromoterActivity(normalizedPromoterCounts.star,
                                                               promoterAnnotationData)
  save(absolutePromoterActivity.star,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivity.star.RData"))
  
  # Calculate gene expression
  geneExpression.star <- getGeneExpression(absolutePromoterActivity.star)
  save(geneExpression.star,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"geneExpression.star.RData"))
  
  # Calculate relative promoter activity
  relativePromoterActivity.star <- getRelativePromoterActivity(absolutePromoterActivity.star,
                                                               geneExpression.star)
  save(relativePromoterActivity.star,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"relativePromoterActivity.star.RData"))
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

################################################
####Merge replicates with mean
################################################
if(TRUE){
  dir.create(file.path(prewd,"ProcessedData","PromoterObjects",method))
  dir.create(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method))
  
  # load AbsolutePromoterActivity and sample info
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivity.star.RData"))
  
  # make the sample info
  sample.info = data.frame(samplename = colnames(absolutePromoterActivity.star)[-1:-2],stringsAsFactors = FALSE)
  sample.info$condition = gsub("_rep.*","", sample.info$samplename)
  # sample.info$tissue = sapply(strsplit(sample.info$samplename, "_"), "[[", 1)
  # sample.info$stage = sapply(strsplit(sample.info$samplename, "_"), "[[", 2)
  
  ##AbsolutePromoterActivity merge replicates with mean 
  absolute.promoter.activity.star.mean = cbind()
  for(sample in unique(sample.info$condition)){
    sample.rownames.selected = sample.info[sample.info$condition==sample,"samplename"]
    absolute.promoter.activity.star.mean = cbind(absolute.promoter.activity.star.mean, rowMeans(absolutePromoterActivity.star[,sample.rownames.selected,drop=FALSE], na.rm = TRUE))
  }
  colnames(absolute.promoter.activity.star.mean) = unique(sample.info$condition)
  save(absolute.promoter.activity.star.mean,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolute.promoter.activity.star.mean.RData"))
  
  # load RelativePromoterActivity
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"relativePromoterActivity.star.RData"))
  ##merge replicates with mean 
  relative.promoter.activity.star.mean = cbind()
  for(sample in unique(sample.info$condition)){
    sample.rownames.selected = sample.info[sample.info$condition==sample,"samplename"]
    relative.promoter.activity.star.mean = cbind(relative.promoter.activity.star.mean, rowMeans(relativePromoterActivity.star[,sample.rownames.selected,drop=FALSE], na.rm = TRUE))
  }
  colnames(relative.promoter.activity.star.mean) = unique(sample.info$condition)
  save(relative.promoter.activity.star.mean,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"relative.promoter.activity.star.mean.RData"))
  
  # load GeneExpression
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"geneExpression.star.RData"))
  ##merge replicates with mean 
  gene.expression.star.mean = cbind()
  for(sample in unique(sample.info$condition)){
    sample.rownames.selected = sample.info[sample.info$condition==sample,"samplename"]
    gene.expression.star.mean = cbind(gene.expression.star.mean, rowMeans(geneExpression.star[,sample.rownames.selected,drop=FALSE], na.rm = TRUE))
  }
  colnames(gene.expression.star.mean) = unique(sample.info$condition)
  save(gene.expression.star.mean,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"gene.expression.star.mean.RData"))
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}  

#######################################################################
#####Major, minor and inactive Promoter classification (STAR alignment)
#######################################################################

if(TRUE){
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivity.star.RData"))
  # absolutePromoterActivity.star = absolutePromoterActivity.star[,1:10]
  absolutePromoterActivity.star[is.na(absolutePromoterActivity.star)] = 0
  absolutePromoterActivity.star$AverageActivity = apply(absolutePromoterActivity.star[,-1:-2],1,mean)
  
  
  major.minor.promoter.classification = function(x){
    gene.id = absolutePromoterActivity.star[,2]
    gene.split = split(x,gene.id)
    
    major.minor.value.change.per.gene = function(z){
      z=as.numeric(z)
      if (max(z)<0.25){
        z[z<0.25]=0 
      } else {
        z[z<0.25]=0 
        z[z==max(z)] =100000
        z[z>=0.25 & z < max(z)]=1000
      }
      z[z==100000] = "Major"
      z[z==1000] = "Minor"
      z[z==0] = "Inactive"
      return(z)
    }
    
    category = lapply(gene.split, major.minor.value.change.per.gene)
    return(as.character(melt(category)[,1])) #convert list to vector without changing order
    
    }
  
  ##classify different types
  absolutePromoterActivity.star$PromoterCategory = as.character(apply(absolutePromoterActivity.star[,"AverageActivity",drop=FALSE],2,major.minor.promoter.classification))
  # absolutePromoterActivity.star.category = cbind(absolutePromoterActivity.star[,1:2], absolutePromoterActivity.star.category)
  # rownames(absolutePromoterActivity.star.category) = rownames(absolutePromoterActivity.star)
  absolutePromoterActivityCategory.star = absolutePromoterActivity.star[,c("promoterId","geneId","AverageActivity","PromoterCategory")]
  save(absolutePromoterActivityCategory.star,file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityCategory.star.RData"))

  rm(list=setdiff(ls(),basic.objects))
  gc()
}

#######################################################################
#####Plot Major, minor and inactive Promoter Percentage (STAR alignment)
#######################################################################

if (TRUE) {
  dir.create(file.path(prewd,"Output","PromoterAnnotation"))
  dir.create(file.path(prewd,"Output","PromoterAnnotation",method))
  dir.create(file.path(prewd,"Output","PromoterAnnotation",method,norm.method))
  
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityCategory.star.RData"))
  absolutePromoterActivityCategory.star.PromoterCategory=as.data.frame(table(absolutePromoterActivityCategory.star$PromoterCategory))
  
  ##Percentage of three promoter category (Major, minor, inactive) 
  absolutePromoterActivityCategory.star.PromoterCategory$Var1 <- factor(absolutePromoterActivityCategory.star.PromoterCategory$Var1, levels = c("Major","Minor","Inactive"))
  # absolutePromoterActivityCategory.star.PromoterCategory$Var2 = "Category"
  colnames(absolutePromoterActivityCategory.star.PromoterCategory)[1] = "Type"
  
  p1<-ggplot(absolutePromoterActivityCategory.star.PromoterCategory, aes(x=1,y=Freq,fill=Type, width=0.5)) +     
    geom_bar(position="stack",stat="identity") + 
    labs(y = "Numbers of promoters") + 
    # scale_x_discrete(c("Inactive","Minor","Major")) +
    scale_fill_manual(values =c("#ff1e56", "#ffac41", "#323232")) +
    scale_color_manual(values =c("#ff1e56", "#ffac41", "#323232")) +
    theme_classic.add +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "right")
  p1
  ggsave(p1, file=file.path(prewd,"Output","PromoterAnnotation",method,norm.method,"promoter.activity.category.percentage.pdf"),width=2.7,height = 2.9,units = "in",dpi=300) #Fig. 1c
  
  ##Comparison of Promoter Activity between three types
  # absolutePromoterActivityCategory.star.plot=absolutePromoterActivityCategory.star[,3:4]
  absolutePromoterActivityCategory.star.minor = absolutePromoterActivityCategory.star[absolutePromoterActivityCategory.star$PromoterCategory %in% c("Minor"),]
  absolutePromoterActivityCategory.star.unique.plot = absolutePromoterActivityCategory.star[absolutePromoterActivityCategory.star$geneId %in% unique(absolutePromoterActivityCategory.star.minor$geneId), ]
  
  head(absolutePromoterActivityCategory.star.unique.plot)  
  absolutePromoterActivityCategory.star.unique.plot$PromoterCategory <- factor(absolutePromoterActivityCategory.star.unique.plot$PromoterCategory, levels = c("Major","Minor","Inactive"))
  p2<-ggplot(absolutePromoterActivityCategory.star.unique.plot, aes(x=PromoterCategory, y=AverageActivity, color=PromoterCategory,fill=PromoterCategory)) + 
    geom_boxplot(alpha=1,outlier.shape = NA, notch =TRUE, width=0.5) +
    scale_fill_manual(values =c("#ff1e56", "#ffac41", "#323232")) + 
    scale_color_manual(values =rep("black",3))+
    labs(x="", y = "Average promoter activity")+
    stat_compare_means(label = "p.signif",method = "wilcox.test",ref.group = "Major",na.rm = TRUE, label.y = max(absolutePromoterActivityCategory.star.unique.plot$AverageActivity) - 1) + 
    ggtitle("")+theme_classic.add+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
  p2
  ggsave(p2,file=file.path(prewd,"Output","PromoterAnnotation",method,norm.method,"promoter.activity.category.comparison.pdf"),width=2.5,height = 2.5,units = "in",dpi=300) #Fig. 1b
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

#######################################################################
#####Plot Major, minor Percentage in different Promoter Positions  (STAR alignment)
#######################################################################

if(TRUE){
  
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityCategory.star.RData"))
  absolutePromoterActivityCategory.star.minor = absolutePromoterActivityCategory.star[absolutePromoterActivityCategory.star$PromoterCategory %in% c("Minor"),]
  absolutePromoterActivityCategory.star.unique.plot = absolutePromoterActivityCategory.star[absolutePromoterActivityCategory.star$geneId %in% unique(absolutePromoterActivityCategory.star.minor$geneId), ]
  absolutePromoterActivityCategory.star.gene.split = split(absolutePromoterActivityCategory.star$PromoterCategory,absolutePromoterActivityCategory.star$geneId)

  ##extract each position for each gene
  active.promoter.position = function(z,number){
    if (any(z %in% "Major")){
      position = z[number]
    } else {position = "None"}
    return(position)
  }
  
  absolutePromoterActivityCategory.star.gene.split.position1 = table(sapply(absolutePromoterActivityCategory.star.gene.split,function(x){active.promoter.position(x,1)}))

  absolutePromoterActivityCategory.star.gene.split.position2 = table(sapply(absolutePromoterActivityCategory.star.gene.split,function(x){active.promoter.position(x,2)}))

  absolutePromoterActivityCategory.star.gene.split.position3 = table(sapply(absolutePromoterActivityCategory.star.gene.split,function(x){active.promoter.position(x,3)}))

  absolutePromoterActivityCategory.star.gene.split.position4 = table(sapply(absolutePromoterActivityCategory.star.gene.split,function(x){active.promoter.position(x,4)}))

  absolutePromoterActivityCategory.star.gene.split.position5 = sapply(absolutePromoterActivityCategory.star.gene.split,function(x){active.promoter.position(x,5:100)})
  absolutePromoterActivityCategory.star.gene.split.position5 = table(unlist(absolutePromoterActivityCategory.star.gene.split.position5))
  
  absolutePromoterActivityCategory.star.gene.split.position.plot=cbind(absolutePromoterActivityCategory.star.gene.split.position1,absolutePromoterActivityCategory.star.gene.split.position2,
                                                                  absolutePromoterActivityCategory.star.gene.split.position3, absolutePromoterActivityCategory.star.gene.split.position4,
                                                                  absolutePromoterActivityCategory.star.gene.split.position5)
  colnames(absolutePromoterActivityCategory.star.gene.split.position.plot)=c(1:5)
  
  absolutePromoterActivityCategory.star.gene.split.position.forplot = as.data.frame(absolutePromoterActivityCategory.star.gene.split.position.plot)[c("Major","Minor"), ]
  absolutePromoterActivityCategory.star.gene.split.position.forplot1 = melt(t(absolutePromoterActivityCategory.star.gene.split.position.forplot))
  plot.width = apply(absolutePromoterActivityCategory.star.gene.split.position.forplot,2,sum)/sum(absolutePromoterActivityCategory.star.gene.split.position.forplot)
  # Stacked
  p1<-ggplot(absolutePromoterActivityCategory.star.gene.split.position.forplot1, aes(x=Var1, y=value, fill=Var2))+ 
    geom_bar(position="fill",stat="identity",width=log10(plot.width*100),colour=rep("black",10))+
    scale_y_continuous(breaks = seq(0,1,by=0.2)) +
    scale_fill_manual(values =c("#ff1e56", "#ffac41")) +
    # scale_color_manual(values =rep("black",10)) +
    theme_classic.add + 
    labs(y = "Proportion of promoter types",x="Promoter position (5' to 3')")+
    theme(legend.title = element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p1
  ggsave(p1,file=file.path(prewd,"Output","PromoterAnnotation",method,norm.method,"promoter.activity.position.category.pdf"),width=3,height = 3,units = "in",dpi=300) #Fig. 1d
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc() 
}

#######################################################################
#####Plot Active Promoter Percentage (STAR alignment)
#######################################################################

if(TRUE){
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"geneExpression.star.RData"))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityCategory.star.RData"))
  
  ##Correlation between Average Activity and Average Expression
  geneExpression.star$AverageExpression<- apply(geneExpression.star[,-1:-2],1,mean)
  absolutePromoterActivityCategory.star.Major = absolutePromoterActivityCategory.star[absolutePromoterActivityCategory.star$PromoterCategory %in% c("Major"),]
  absolutePromoterActivityCategory.star.Major.merge.geneExpression = merge(absolutePromoterActivityCategory.star.Major,geneExpression.star[,c("geneId","AverageExpression")],by.x="geneId",by.y="geneId")
  
  p1=ggplot(absolutePromoterActivityCategory.star.Major.merge.geneExpression, aes(x=AverageExpression, y=AverageActivity)) + 
    geom_point(size=0.05,colour="black",alpha=0.9)+
    xlab("Average gene expression") +
    ylab("Average promoter activity") +
    theme_classic.add +xlim(0,15)+ylim(0,15)
  p1  
  ggsave(p1,file=file.path(prewd,"Output","PromoterAnnotation",method,norm.method,"promoter.activity.geneexpression.correlation.pdf"),width=3,height = 3,units = "in",dpi=300)
  
  
  ##single or multiple promoter per gene
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityCategory.star.RData"))
  absolutePromoterActivityCategory.star.gene.split = split(absolutePromoterActivityCategory.star$PromoterCategory,absolutePromoterActivityCategory.star$geneId)
  
  ##Percentage of three promoter category (Major, minor, inactive) plot
  absolutePromoterActivityCategory.star.PromoterCategory=as.data.frame(table(absolutePromoterActivityCategory.star$PromoterCategory))
  absolutePromoterActivityCategory.star.PromoterCategory$Var1 <- factor(absolutePromoterActivityCategory.star.PromoterCategory$Var1, levels = c("Major","Minor","Inactive"))
  absolutePromoterActivityCategory.star.PromoterCategory$Freq <- c(length(absolutePromoterActivityCategory.star.gene.split)-length(unlist(sapply(absolutePromoterActivityCategory.star.gene.split,function(x){return(grep("Major",x))})))-length(unlist(sapply(absolutePromoterActivityCategory.star.gene.split,function(x){return(grep("Minor",x))}))),
                                                                   length(unlist(sapply(absolutePromoterActivityCategory.star.gene.split,function(x){return(grep("Major",x))}))),
                                                                   length(unlist(sapply(absolutePromoterActivityCategory.star.gene.split,function(x){return(grep("Minor",x))})))
                                                                   )
  colnames(absolutePromoterActivityCategory.star.PromoterCategory)[1] = "Type"
  
  p1<-ggplot(absolutePromoterActivityCategory.star.PromoterCategory[c(2,3,1),], aes(x=1,y=Freq,fill=Type, width=0.5)) +     
    geom_bar(position="stack",stat="identity") + 
    labs(y = "Numbers of promoters") + 
    # scale_x_discrete(limits = c("Major","Minor","Inactive")) +
    scale_fill_manual(values =c("#ff1e56", "#ffac41", "#323232")) +
    scale_color_manual(values =c("#ff1e56", "#ffac41", "#323232")) +
    theme_classic.add +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "right")
  p1
  ggsave(p1,file=file.path(prewd,"Output","PromoterAnnotation",method,norm.method,"promoter.activity.category.percentage.genewise.pdf"),width=2.7,height = 3,units = "in",dpi=300)
  
  
  ##Percentage of active promoter distribution plot
  active.promoter.classification = function(z){
    if (any(z %in% "Major")){
      
      if(any(z%in% "Minor")){
        category = "Multipromoter.Multiactive"
      } else if ((!any(z %in% "Minor")) & any(z %in% "Inactive")){
        category = "Multipromoter.Singleactive"
      } else if ((!any(z %in% "Minor")) & (!any(z %in% "Inactive"))){
        category = "Singlepromoter.Singleactive"
      }
      
    } else {category = "Inactive"}
    
    return(category)
  }
  
  absolutePromoterActivityCategory.star.gene.split.category = sapply(absolutePromoterActivityCategory.star.gene.split,active.promoter.classification)
  
  absolutePromoterActivityCategory.star.gene.split.category.plot=as.data.frame(table(absolutePromoterActivityCategory.star.gene.split.category))
  colnames(absolutePromoterActivityCategory.star.gene.split.category.plot)=c("Category","Freq")
  absolutePromoterActivityCategory.star.gene.split.category.plot$Category=c("Inactive","Multiple active promoters\n(Multi promoter genes)","Single active promoters\n(Multi promoter genes)","Single active promoters\n(Single promoter genes)")
  
  p1<-ggplot(absolutePromoterActivityCategory.star.gene.split.category.plot[-1, ], aes(x=1,y= Freq/sum(Freq)*100,fill=Category)) + 
    geom_bar(position="stack",stat="identity") + 
    labs(y = "% of active promoters")+ coord_flip()+  
    scale_fill_manual(values =c("#ff6464", "#db3056", "#851d41")) + 
    scale_color_manual(values =rep("black",3))+
    theme_classic.add+ 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
          legend.position = "top",legend.direction="vertical",
          legend.title = element_blank()) #,text = element_text(size = 24)  
  p1
  ggsave(p1,file=file.path(prewd,"Output","PromoterAnnotation",method,norm.method,"promoter.activety.single.multiple.category.pdf"),width=2.5,height =2.5,units = "in",dpi=300) #Fig. 1f
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc() 
}


################################################
####Plotting promoter annotations
################################################
if(TRUE){
  load(file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
  promoterCoordinates(promoterAnnotationData)[is.na(promoterCoordinates(promoterAnnotationData)$internalPromoter)]$internalPromoter = FALSE
  promoterCoordinates(promoterAnnotationData) = promoterCoordinates(promoterAnnotationData)[!promoterCoordinates(promoterAnnotationData)$internalPromoter]
  
  promoter.id.mapping.data<-promoterIdMapping(promoterAnnotationData)
  promoter.id.mapping.data.geneid.table<-as.data.frame(table(unique(promoter.id.mapping.data[,c("promoterId","geneId")])$geneId))
  promoter.id.mapping.data.geneid.table[promoter.id.mapping.data.geneid.table$Freq>11,"Freq"]=11
  promoter.id.mapping.data.geneid.table.plot = as.data.frame(table(promoter.id.mapping.data.geneid.table$Freq))
  
  dir.create(file.path(prewd,"Output"))
  dir.create(file.path(prewd,"Output","PromoterAnnotation"))
  
  ##Promoter number hist distribution
  p1<-ggplot(promoter.id.mapping.data.geneid.table.plot, aes(x=Var1,y=Freq)) + 
    geom_bar(stat="identity",color="black", fill="white") + 
    geom_text(aes(label=Freq), vjust=-0.2, size=2.5)+ 
    scale_x_discrete(labels=c(1:10,">=11")) +
    labs(x="Promoter number", y = "Number of genes")+
    ggtitle("Promoter number distribution")+theme_classic.add+theme(plot.title = element_text(hjust = 0.5))
  p1
  ggsave(p1,file=file.path(prewd,"Output","PromoterAnnotation",'promoter.number.hist.all.pdf'),width=3.5,height =3,units = "in",dpi=300) ##Fig. S1a
  
  ##Promoter number hist distribution without single promoter
  promoter.id.mapping.data.geneid.table.plot2 = promoter.id.mapping.data.geneid.table.plot[-1,]
  p2<-ggplot(promoter.id.mapping.data.geneid.table.plot2, aes(x=Var1,y=Freq)) + 
    geom_bar(stat="identity",color="black", fill="white") + 
    geom_text(aes(label=Freq), vjust=-0.2, size=2.5)+
    scale_x_discrete(labels=c(2:10,">=11")) +
    labs(x="Promoter number", y = "Number of genes")+
    ggtitle("Promoter number distribution")+theme_classic()+theme(plot.title = element_text(hjust = 0.5))
  p2
  ggsave(p2,file=file.path(prewd,"Output","PromoterAnnotation",'promoter.number.hist.without1.pdf'),width=3.5,height =3,units = "in",dpi=300)
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}




