prewd="Mammalian2019Promoter/CM2019RNAseqMouse"
stage.order = c("E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5", "E17.5", "E18.5", "P0", "P03", "P14", "P28", "P63")
sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects","sample.sheet.txt"),header=TRUE, row.names = 1, stringsAsFactors = FALSE)
rownames(sample.info) = sample.info$newnamegender

library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(GenomicFeatures)
library(ensembldb)
library(GenomicAlignments)
library(AnnotationDbi)
library(dplyr)
#library(proActiv)
library(ggplot2)
library(reshape2)
#library(maSigPro)
library(doParallel)
library(foreach)
library(pheatmap)
#library(webshot)
library(plotly)

method="Proactiv"
# method="Dexseq"
norm.method = "edger"
# norm.method = "deseq2"
gender = "Male"
# gender = "Female"
# gender = "ale"

source("~/scripts/tjcm.promoter.stage.specific.function.R")

if(gender == "Male"){
  tissue.order = c("Brain", "Cerebellum", "Heart", "Kidney" , "Liver", "Testis")
} else if (gender == "Female") {
  tissue.order = c("Ovary") #"Brain", "Cerebellum", "Heart", "Kidney" , "Liver", 
} else {
  tissue.order = c("Brain", "Cerebellum", "Heart", "Kidney" , "Liver", "Ovary", "Testis")
}
mc.cores= 20
theme_classic.add = theme_classic()+ theme(legend.position = "none",strip.background = element_blank(),strip.text = element_blank(), plot.title = element_text(hjust = 0.5, family = "Helvetica"),axis.text = element_text(colour = "black", family = "Helvetica"), axis.title=element_text(colour = "black", family = "Helvetica"), text = element_text(colour = "black", family = "Helvetica"))
theme_bw.add = theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, family = "Helvetica"),axis.text = element_text(colour = "black", family = "Helvetica"), axis.title=element_text(colour = "black", family = "Helvetica"), text = element_text(colour = "black", family = "Helvetica"))
basic.objects=c("basic.objects",ls(),"name")

######################################################
##regression fit to calculate pvalue and R2 correlation for each promoter 
##IsoModel 1. only consider genes with more than one expressed(pre-filtered sum> 10, any>1) isoforms, calculate p-value to find differentially expressed genes (I skipped this pvalue filteration)
##2. p.vector: fit a linear model and calculate p-value and p-adjusted to find significant isoforms
##3. Tfit: calculate R correlation to filter significant genes, and test hypothesis to find group difference on top of time difference (if there are several groups), 
##   here is to find differences between experimental groups, whether the treatment is significantly different from the control
##Ouput 1. data: all the original data
##2. gen: all the isoform corresponding gene name
##3. design: design matrix
##4. DSG: differentially expressed genes (if we do not filter, then it is all genes with more than 1 promoter)
##5. pvector result: fit$G,g # returns the number of all isoforms
##                   fit$dat # is a matrix with the all isoforms and their expression values  
##                   fit$p.vector # gives p-value at the Q false discovery control level for all isoforms (fit the linear model)
##                   fit$p.adjusted # gives p-adjusted at the Q false discovery control level for all isoforms
##                   fit$i # returns the number of significant isoforms (p-adjusted <= Q)
##                   fit$SELEC = fit$dat[which(fit$p.adjusted<=0.2),] is a matrix with the significant isoforms and their expression values (p-adjusted <= Q) 
##                   fit$dis # the design matrix for time
##                   fit$groups.vector # group info for the design matrix columns
##                   fit$edesign # the design matrix for time
##                   fit$family # the model being used, NB or Gaussian
##                   fit$min.obs # genes with less than this number of true numerical values will be excluded from the analysis. Minimum value to estimate the model is (degree+1)xGroups+1. Default is 6.
##                   fit$Q # the p-adjusted cutoff that you set
##6. Tfit result:    tfit$G # returns the number of all isoforms
##                   tfit$g # returns the number of significant isoforms from pvector result
##                   tfit$dat # is a matrix with the significant isoforms from pvector result and their expression values  
##                   sol: table 3 in first paper; p-value of the regression coefficients of the selected variables (time, time2, coldvsControl, saltvsControl); first p-value is the model compared to the mean; 
##                   sol: beta0,1,2 means beta0, sigma0, delta0, which are the regression p-value corresponding to the reference group. bi, di, gi.. .li are the regression coefficients that account for specific differences (linear, quadratic, cubic, etc.) between the (i + 1)-th group profile and the first group (reference) profil 
##                   coefficients: table 3 in first paper; beta0,Time,Time2 means beta0, sigma0, delta0, coefficients (slope) of the regression coefficients of the selected variables (time, time2, coldvsControl, saltvsControl);
##                   t.score: table 3 in first paper; t.score of the regression coefficients of the selected variables (time, time2, coldvsControl, saltvsControl);
##                   sig.profiles is a matrix with the significant isoforms and their expression values (p-adjusted <= Q)
##                   group.coeffs: calculate beta0 and beta coefficent for different groups; Coldbeta0 = beta0+betaColdvsControl; Coldbeta1 = betaTime+betaTimexCold; when there is only 1 group, it is the same as coefficients
##                   variables: time variable in design matrix
##                   dis # the design matrix for time
##                   groups.vector # group info for the design matrix columns
##                   edesign # the design matrix for time
##                   step.method: forward or backward 
##get.siggenes:      select only significantly changed genes. always use vars="groups", significant.intercept = "dummy", because it exclude beta0 (y-intercept) difference, only keeps difference in terms of time
##                   vars = c("all", "each", "groups"); all (beta0, time, time2) means any p-value <= 0.05 is fine; each means seperate each variable; group means seperate each condition(salt, cold)
##                   group: significant.intercept = c("all", "dummy", "none"); all means anyone is not NA ("p.valor_beta0" "p.valor_Time"  "p.valor_Time2"); dummy means any two is not NA ("p.valor_Time"  "p.valor_Time2"); none is same as dummy if you only have time variable
##Output: get significantly changed promoters
##1. Model:          all isomodel (p.vector, Tfit) result from the last step
##2. get2            get.siggenes result
##                   get2$sig.genes$group ##selected signicantly isoforms, all Tfit result
##                   get2$summary # names of selected signicantly isoforms
##4. DET:            differentially expressed isoforms; result from get2$summary
##3. DSG:            differentially expressed genes correspond to differentially expressed isoforms
##3. List0:          genes that are not correspond to differentially expressed isoforms from last step, since we did R correlation filter in this step.
##4. NumIso.by.gene: how many differentially expressed isoforms in one gene, we are interested into those genes with more than 1 changing isoforms.
##Seegenes:          clustering DE isoforms into clusters
##                   cluster.all: take only significantly changed genes 
##Output:            Model: the isomodel (p.vector, Tfit) result which it used as an input
##                   get2: get.siggenes result which it used as an input
##                   NumIso.by.gene: how many differentially expressed isoforms in one gene, we are interested into those genes with more than 1 changing isoforms.
##                   cut: cluster for each isoform
##                   names.genes: each isoform corresponding gene name 
##TableDS:           excluding genes with only 1 isoform changing; make an interaction table about two/more isoforms in a gene 
##Output:            IsoMajorMinor: major and minor promoter of the gene
##                   IsoClusters: Clusters of major and minor promoter
##                   IsoTable: isoforms interaction table
##PodiumChange:      major isoforms changes; select genes with major isoforms changed in any time points or genes with major isoforms changed between two groups in any time points. 
##                   only.sig.iso: take only significantly changed isoforms from get.siggenes, otherwise take correspond genes of all isoforms from get.siggenes, then use all isoforms of this gene
##                   comparison = c("any", "groups","specific"); "any": select genes with major isoforms changed in any time points in any group; "groups": genes with major isoforms changed between two groups in any time points; "specific": genes with major isoforms changed between two specific groups at a specific time point
##Output:            L: selected uniquely genes with major isoforms changed at any time point
##                   data.L: genes corresponding isoforms expression level
##                   gen.L: isoforms corresponding gene name
############Two important table afterwards#####################
##cluster category (up,down,flat) info: significant.isoform.tissue.see.cluster,cluster0
##table category info: significant.isoform.tissue.table.category ("Down_Down","Down_Flat","Flat_Down","Down_Up","Up_Down","Flat_Up","Up_Flat","Up_Up")
##              


#################################################################################
###Time series analysis for different developmental stages in each tissue
#################################################################################
if(TRUE){ #file.exists(file.path(prewd,"ProcessedData","PromoterObjects","human.sample.sheet.txt"))
  ###different promoters contributed difference along the development
  dir.create(file.path(prewd,"ProcessedData","StageObjects"))
  dir.create(file.path(prewd,"ProcessedData","StageObjects",method))
  dir.create(file.path(prewd,"ProcessedData","StageObjects",method,norm.method))
  dir.create(file.path(prewd,"Output","StageSpecific"))
  dir.create(file.path(prewd,"Output","StageSpecific",method))
  dir.create(file.path(prewd,"Output","StageSpecific",method,norm.method))
  
  ####################################
  ## load normalized counts information
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivity.star.RData"))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"normalizedPromoterCounts.star.RData"))
  absolutePromoterActivity.star = cbind(absolutePromoterActivity.star[,1:2],normalizedPromoterCounts.star[rownames(absolutePromoterActivity.star),])
  
  ## select gender
  absolutePromoterActivity.star = absolutePromoterActivity.star[,c("geneId", grep(gender, colnames(absolutePromoterActivity.star), value=TRUE) ), drop=FALSE]
  
  ## select highly expressed ones
  absolutePromoterActivity.combine = absolutePromoterActivity.star[rowSums(absolutePromoterActivity.star[,-1])>=50 & 
                                                                              apply(absolutePromoterActivity.star[,-1],1,function(x){max(x)>1}),]
  
  ####################################
  ##load sample info 
  # sample.info = data.frame(samplename = colnames(absolutePromoterActivity.star[,-1]), stage = unlist(sapply(strsplit(colnames(absolutePromoterActivity.star[,-1]),"_"), function(x){x[2]})), tissue = gsub("_.*","",colnames(absolutePromoterActivity.star[,-1])), condition = unlist(sapply(strsplit(colnames(absolutePromoterActivity.star[,-1]),"_"), function(x){paste(x[1],x[2], sep = "_")})), stringsAsFactors = FALSE)
  sample.info.combine = sample.info[grep(gender, sample.info$gender),]
  
  ##time
  if(length(grep("p",sample.info.combine$stage)) != 0){
    time.change = gsub("6month","1.5",gsub("0month", "0.76", gsub("year","",gsub("week", "", sample.info.combine$stage))))  #1:length(unique(sample.info.combine$stage))
    time = 1:length(time.change)
    time[grep("b",time.change)] = log2(as.numeric(gsub("b", "", time.change[grep("b",time.change)])) * 7)
    time[grep("p",time.change)] = log2(as.numeric(gsub("p", "", time.change[grep("p",time.change)])) * 365)
  } else if (length(grep("P",sample.info.combine$stage)) != 0){
    time.change = sample.info.combine$stage  #1:length(unique(sample.info.combine$stage))
    time = 1:length(time.change)
    time[grep("E",time.change)] = log2(as.numeric(gsub("E", "", time.change[grep("E",time.change)])))
    time[grep("P",time.change)] = log2(as.numeric(gsub("P", "", time.change[grep("P",time.change)])) + 21)
  }
  sample.info.combine$time = time #[sample.info.combine$stage]
  sample.info.combine = sample.info.combine[order(sample.info.combine$tissue,sample.info.combine$time),]
  
  ##replicate
  replicate = 1:length(unique(sample.info.combine$condition)) #paste(sample.info.combine$stage,sample.info.combine$tissue,sep="_")
  names(replicate) = unique(sample.info.combine$condition)
  sample.info.combine$replicatenew = replicate[sample.info.combine$condition]  #unlist(sapply(split(sample.info.combine$stage,factor(replicate,levels = unique(replicate))),function(x)return(1:length(x))))
  sample.info.combine.sup = as.data.frame(model.matrix(~ 0 + tissue, sample.info.combine))
  colnames(sample.info.combine.sup) = gsub("tissue","",colnames(sample.info.combine.sup))
  sample.info.all = cbind(sample.info.combine,sample.info.combine.sup)
  sample.info.combine = sample.info.all[,!colnames(sample.info.all) %in% colnames(sample.info)]
  colnames(sample.info.combine) = gsub("new","",colnames(sample.info.combine))
  
}  
  
##time differences comparing one tissue to others 
########################################################
for (name in tissue.order){ #"Brain", "Cerebellum", "Heart", "Ovary", "Kidney", "Liver" c("Testis")
  message(name)
  dir.create(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name))
  dir.create(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender))
  dir.create(file.path(prewd,"Output","StageSpecific",method,norm.method,name))
  dir.create(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender))

  sample.info.combine=sample.info.combine[,c("time","replicate",name, grep(name,tissue.order,invert = TRUE, value = TRUE))]
  design.all <- make.design.matrix(sample.info.combine, degree = 3)
  isoform.model.all <- IsoModel(data=absolutePromoterActivity.combine[,-1], gen=absolutePromoterActivity.combine[,1], design=design.all, Q=1, Qfit = 0.05, counts=TRUE, theta = 10, step.method = "backward") #"forward", Q = 0.05, min.obs=15, minorFoldfilter = 10,
  save(isoform.model.all, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("isoform.model.all.RData")))
}

##time differences within one tissue 
####################################
for (name in tissue.order){ #"Brain", "Cerebellum", "Heart", "Ovary", "Kidney", "Liver" c("Testis")
  message(name)
  sample.info.tissue.combine=sample.info.combine[grep(name, rownames(sample.info.combine)),c("time","replicate",name)]
  sample.info.tissue.combine$replicate = sample.info.tissue.combine$replicate - min(sample.info.tissue.combine$replicate) + 1
  colnames(sample.info.tissue.combine) = c("time", "replicate", "group")
  
  ####################################
  ##Counts selected for one tissue
  absolutePromoterActivity.tissue = absolutePromoterActivity.star[,c("geneId",rownames(sample.info.tissue.combine))]
  absolutePromoterActivity.tissue.combine = absolutePromoterActivity.tissue[rowSums(absolutePromoterActivity.tissue[,-1])>=10 &
                                                                              apply(absolutePromoterActivity.tissue[,-1], 1, function(x){max(x)>1}),]

  design.tissue <- make.design.matrix(sample.info.tissue.combine, degree = 3)
  isoform.model.tissue <- IsoModel(data=absolutePromoterActivity.tissue.combine[,-1], gen=absolutePromoterActivity.tissue.combine[,1], design=design.tissue, Q=1, Qfit = 0.05, counts=TRUE, theta = 10, step.method = "backward") #"forward", Q = 0.05, min.obs=15, minorFoldfilter = 10, 
  save(isoform.model.tissue, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("isoform.model.tissue.RData")))
  # isoform.model.tissue.all <- IsoModel(data=absolutePromoterActivity.tissue.combine[,-1], gen=absolutePromoterActivity.tissue.combine[,1], design=design.tissue, Q=1, Qfit = 1, counts=TRUE, theta = 10, step.method = "backward") #, Q = 0.05, min.obs=15, minorFoldfilter = 10,
  # save(isoform.model.tissue.all, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("isoform.model.tissue.all.RData")))
}

rm(list=setdiff(ls(),basic.objects))
gc()

##significantly changing isoforms analysis per tissue 
####################################################
for (name in tissue.order){ #"Brain", "Cerebellum", "Heart", "Ovary", "Kidney", "Liver" c("Testis")
  message(name)
  ##significantly changing isoforms in terms of time per tissue 
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("isoform.model.tissue.RData")))
  significant.isoform.tissue <- getDS(Model = isoform.model.tissue, rsq = 0.4, vars = "all") #groups, all, each, we should always use groups, but I checked the result in this dataset, the result of all is the same as groups
  length(significant.isoform.tissue$List0)
  length(significant.isoform.tissue$DSG)
  length(significant.isoform.tissue$DET)
  save(significant.isoform.tissue, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
  
  ##clustering DE isoforms to see if developmentally down or up regulated
  pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("seegenes.cluster.pdf")))
  print(significant.isoform.tissue.see <- seeDS(get = significant.isoform.tissue, cluster.all=FALSE, k=2) ) #cluster.all here refers to genes with more than 1 isoforms
  dev.off()
  save(significant.isoform.tissue.see, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.RData")))
  significant.isoform.tissue.see.cluster = data.frame(cluster=significant.isoform.tissue.see$cut,geneid=significant.isoform.tissue.see$names.genes, numiso = significant.isoform.tissue.see$NumIso.by.gene[significant.isoform.tissue.see$names.genes], stringsAsFactors = FALSE)
  save(significant.isoform.tissue.see.cluster, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster.RData")))
  write.table(significant.isoform.tissue.see.cluster, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster.txt")), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
  
  ##In one gene, only one isoform is developmental dynamically changing while the others are not, extract the others, naming cluster0
  promoter.gene.correspond = data.frame(promoterid = rownames(significant.isoform.tissue$Model$data), geneid = significant.isoform.tissue$Model$gen, stringsAsFactors = FALSE)
  significant.isoform.tissue.transcript.number.singlechange = significant.isoform.tissue.see.cluster[significant.isoform.tissue.see.cluster$numiso==1,"geneid"]
  significant.isoform.tissue.transcript.number.singlechange.all.isoforms = promoter.gene.correspond[promoter.gene.correspond$geneid %in% significant.isoform.tissue.transcript.number.singlechange,]
  significant.isoform.tissue.transcript.number.singlechange.all.isoforms$cluster = 0
  significant.isoform.tissue.see.cluster0 = merge(significant.isoform.tissue.transcript.number.singlechange.all.isoforms, 
                                                  as.data.frame(table(significant.isoform.tissue.transcript.number.singlechange.all.isoforms$geneid)), by.x="geneid", by.y="Var1")
  rownames(significant.isoform.tissue.see.cluster0) = significant.isoform.tissue.see.cluster0$promoterid
  significant.isoform.tissue.see.cluster0 = significant.isoform.tissue.see.cluster0[!rownames(significant.isoform.tissue.see.cluster0) %in% rownames(significant.isoform.tissue.see.cluster),c("cluster","geneid","Freq")]
  colnames(significant.isoform.tissue.see.cluster0)[3] = "numiso"
  save(significant.isoform.tissue.see.cluster0, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster0.RData")))
  
  ##put cluster0 also into significant.isoform.tissue.see
  significant.isoform.tissue.see0 = significant.isoform.tissue.see
  significant.isoform.tissue.see0$get2$sig.genes$sig.profiles = rbind(significant.isoform.tissue.see0$get2$sig.genes$sig.profiles,
                                                                      significant.isoform.tissue.see0$Model$data[rownames(significant.isoform.tissue.see.cluster0),])
  significant.isoform.tissue.see0$get2$summary = c(significant.isoform.tissue.see0$get2$summary, 
                                                   rownames(significant.isoform.tissue.see.cluster0))
  significant.isoform.tissue.see0$cut = c(significant.isoform.tissue.see0$cut, 
                                          setNames(significant.isoform.tissue.see.cluster0$cluster, rownames(significant.isoform.tissue.see.cluster0)) )
  significant.isoform.tissue.see0$names.genes = c(significant.isoform.tissue.see0$names.genes, 
                                                  significant.isoform.tissue.see.cluster0$geneid )
  significant.isoform.tissue.see0$NumIso.by.gene[unique(significant.isoform.tissue.see.cluster0$geneid)] = unique(significant.isoform.tissue.see.cluster0[,c("geneid","numiso")])$numiso
  save(significant.isoform.tissue.see0, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see0.RData")))
  
  ##find isoforms belong to same gene but belong to different clusters
  significant.isoform.tissue.table <- tableDS(seeDS = significant.isoform.tissue.see0)
  print(significant.isoform.tissue.table$IsoTable)
  save(significant.isoform.tissue.table, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.RData")))
  
  ##find genes with major isoform changes at any time point
  significant.isoform.tissue.majoriso.change = PodiumChange(get = significant.isoform.tissue, only.sig.iso=FALSE, comparison="any") #only.sig.iso: take only significantly changed isoforms, otherwise take all isoform correspond genes, then all isoforms of this gene
  length(significant.isoform.tissue.majoriso.change$L)
  save(significant.isoform.tissue.majoriso.change, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.majoriso.change.RData")))
  
  ##significantly changing isoforms comparing this tissue to others
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("isoform.model.all.RData")))
  significant.isoform.all = get.siggenes(tstep = isoform.model.all$Tfit, rsq = 0.1, vars = "groups") #groups, all, each
  save(significant.isoform.all, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.all.RData")))
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

flat.cluster = c(0)
up.cluster = c(2) #2,3,6
down.cluster = c(1) #1,4,5
theme_classic.add = theme_classic()+ theme(legend.position = "none",strip.background = element_blank(),strip.text = element_blank(), plot.title = element_text(hjust = 0.5, family = "Helvetica"),axis.text = element_text(colour = "black", family = "Helvetica"), axis.title=element_text(colour = "black", family = "Helvetica"), text = element_text(colour = "black", family = "Helvetica")) #axis.text = element_text(colour = "black"), axis.title=element_text(colour = "black"), 
basic.objects = c(basic.objects,"up.cluster","down.cluster","flat.cluster","theme_classic.add")

################################################
######isoforms lineplot and heatmap for 5 categories (UpUp, UpDown, DownDown) 
################################################
for (name in tissue.order){#"Brain" c("Testis")
  message(name)
    
  if(TRUE){
    # name = "Testis" "Brain" 
    if(name == "Kidney"){up.cluster = c(1); down.cluster = c(2)} else {up.cluster = c(2); down.cluster = c(1)}
    
    ##load the significant.isoform.tissue and expression info
    # load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.majoriso.change.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see0.RData")))
    load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolute.promoter.activity.star.mean.RData"))  
    absolute.promoter.activity.star.mean = absolute.promoter.activity.star.mean[,grep(gender, colnames(absolute.promoter.activity.star.mean))]
    load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"gene.expression.star.mean.RData"))
    gene.expression.star.mean = gene.expression.star.mean[,grep(gender, colnames(gene.expression.star.mean))]
    
    ##profileplot and heatmap for major and corresponding minor isoforms of different categories
    ###############################################################################################################
    build.cluster.plot = function(cluster1.name="Up_Up", cluster1 = expand.grid(up.cluster, up.cluster)){
      ##gene name belong to this category 
      build.extract.cluster.name = function(cluster){
        y = apply(cluster, 1, function(x){
          getDSPatterns(significant.isoform.tissue.table, x[1], x[2])})
        return(unlist(y))}
      cluster1.genename = build.extract.cluster.name(cluster1)
      # cluster1.genename = cluster1.genename[cluster1.genename %in% significant.isoform.tissue.majoriso.change$L]
      
      ##Major and minor promoter name belong to genes of this category
      major.cluster1 = as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,"Major"])
      build.strsplit = function(x){return(unique(unlist(strsplit(x,"_"))))}
      minor.cluster1 = build.strsplit(as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,"Minor"]))
      write.table(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,], file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.",cluster1.name,".txt")), sep="\t", col.names = NA, row.names = TRUE, quote = FALSE)
      
      ##Expression of Major and minor promter name (mean of replicates)
      significant.isoform.tissue.expression.major.cluster1 = as.data.frame(absolute.promoter.activity.star.mean[c(major.cluster1),grep(name, colnames(absolute.promoter.activity.star.mean))]) #significant.isoform.tissue.see0$get2$sig.genes$sig.profiles[major.cluster1,]
      # tmp.remove = significant.isoform.tissue.expression.major.cluster1[rowSums(significant.isoform.tissue.expression.major.cluster1[,1:11])==0, ] #for tetis only 
      # tmp.remove.table = significant.isoform.tissue.table$IsoMajorMinor[!significant.isoform.tissue.table$IsoMajorMinor$Major %in% rownames(tmp.remove),]
      significant.isoform.tissue.expression.minor.cluster1 = as.data.frame(absolute.promoter.activity.star.mean[c(minor.cluster1),grep(name, colnames(absolute.promoter.activity.star.mean))]) #significant.isoform.tissue.see0$get2$sig.genes$sig.profiles[minor.cluster1,]
      gene.expression.star.mean.tissue = as.data.frame(gene.expression.star.mean[cluster1.genename,grep(name, colnames(gene.expression.star.mean))])
      ##########bootstrap
      build.bootstrap.value = function(input.data){
        size.number = ifelse(length(input.data) <= 100,length(input.data) - 10, 100)
        build.bootstrap = function(input.data,j) {
          random.sample = sample(1:length(input.data),size = size.number, replace = FALSE)
          x.auto.ratio.bootstrap = mean(input.data[random.sample],na.rm = TRUE)
          return(x.auto.ratio.bootstrap)
        }
        bootstrap.value = sapply(1:1000, function(x){build.bootstrap(input.data,x)})
        return(c(max(bootstrap.value), min(bootstrap.value)))
      }
      registerDoParallel(cores = 1)
      significant.isoform.tissue.expression.major.cluster1.bootstrap = foreach(select.col = 1:(ncol(significant.isoform.tissue.expression.major.cluster1)))  %dopar%  build.bootstrap.value(input.data=significant.isoform.tissue.expression.major.cluster1[,select.col])
      significant.isoform.tissue.expression.minor.cluster1.bootstrap = foreach(select.col = 1:(ncol(significant.isoform.tissue.expression.minor.cluster1)))  %dopar%  build.bootstrap.value(input.data=significant.isoform.tissue.expression.minor.cluster1[,select.col])
      gene.expression.star.mean.tissue.bootstrap = foreach(select.col = 1:(ncol(gene.expression.star.mean.tissue)))  %dopar%  build.bootstrap.value(input.data=gene.expression.star.mean.tissue[,select.col])
      ##########combine mean,bootstrap together
      significant.isoform.tissue.expression.major.cluster1.plot = data.frame(Type="Major", mean = colMeans(significant.isoform.tissue.expression.major.cluster1),t(sapply(significant.isoform.tissue.expression.major.cluster1.bootstrap,c)) )
      significant.isoform.tissue.expression.minor.cluster1.plot = data.frame(Type="Minor", mean = colMeans(significant.isoform.tissue.expression.minor.cluster1),t(sapply(significant.isoform.tissue.expression.minor.cluster1.bootstrap,c)) )
      gene.expression.star.mean.tissue.plot = data.frame(Type="Gene", mean = colMeans(gene.expression.star.mean.tissue),t(sapply(gene.expression.star.mean.tissue.bootstrap,c)) )
      significant.isoform.tissue.expression.plot = rbind(significant.isoform.tissue.expression.major.cluster1.plot, significant.isoform.tissue.expression.minor.cluster1.plot, gene.expression.star.mean.tissue.plot)
      colnames(significant.isoform.tissue.expression.plot) = c("Type","mean","max","min")
      significant.isoform.tissue.expression.plot$stage = sapply(strsplit(rownames(significant.isoform.tissue.expression.plot),"_"), "[", 2)
      significant.isoform.tissue.expression.plot$stage = factor(significant.isoform.tissue.expression.plot$stage, levels = stage.order)
      significant.isoform.tissue.expression.plot$Type = factor(significant.isoform.tissue.expression.plot$Type, levels = c("Major","Minor","Gene"))
      
      ##Plot Major and minor expression profile in one metagene lineplot 
      p1 = ggplot(significant.isoform.tissue.expression.plot[significant.isoform.tissue.expression.plot$Type!="Gene",], aes(x=stage, y=mean, color = Type, fill = Type, ymin = min, ymax = max, group = Type)) +
        geom_ribbon(alpha = 0.3, color="white") + 
        geom_point() +
        geom_line() +
        scale_fill_manual(values =c("#ff1e56", "#ffac41")) + 
        scale_color_manual(values =c("#ff1e56", "#ffac41"))+
        ggtitle(paste("Major:",nrow(significant.isoform.tissue.expression.major.cluster1),"; Minor:",nrow(significant.isoform.tissue.expression.minor.cluster1)) ) +
        ylab("") + xlab("") +
        scale_y_continuous(limits = c(0,8.5)) +
        # ylab("log2(Promoter Counts)")+
        theme_classic.add + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) #legend.position="none"
      p1
      ggsave(p1, file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.bootstrap.",cluster1.name,".pdf")), width=2.2, height=2.5)  #Fig. 2a-h, S2-7 middle pannel
      #plus gene
      p2 = ggplot(significant.isoform.tissue.expression.plot, aes(x=stage, y=mean, color = Type, fill = Type, group = Type)) +
        geom_point() +
        geom_line() +
        scale_fill_manual(values =c("#ff1e56", "#ffac41","black")) + 
        scale_color_manual(values =c("#ff1e56", "#ffac41","black"))+
        ggtitle(paste("Major:",nrow(significant.isoform.tissue.expression.major.cluster1),"; Minor:",nrow(significant.isoform.tissue.expression.minor.cluster1)) ) +
        ylab("") + xlab("") +
        scale_y_continuous(limits = c(0,15)) +
        # ylab("log2(Normalized Counts)")+
        theme_classic.add + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) #legend.position="none"
      p2
      ggsave(p2, file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.lineplot.",cluster1.name,".pdf")), width=2.2, height=2.5)  #Fig. 3
      
      ##Plot Major and minor isoform heatmap
      absolute.promoter.activity.mean.tissue.cluster1 = as.data.frame(absolute.promoter.activity.star.mean[c(major.cluster1, minor.cluster1),grep(name, colnames(absolute.promoter.activity.star.mean))])
      colnames(absolute.promoter.activity.mean.tissue.cluster1) = sapply(strsplit(colnames(absolute.promoter.activity.mean.tissue.cluster1),"_"), "[", 2)
      absolute.promoter.activity.mean.tissue.cluster1$category = c(rep("major",length(major.cluster1)), rep("minor",length(minor.cluster1)))
      color.selected= c(major = "#ff1e56", minor = "#ffac41")
      
      pheatmap(absolute.promoter.activity.mean.tissue.cluster1[,-ncol(absolute.promoter.activity.mean.tissue.cluster1)], 
               color = colorRampPalette(c("#77d8d8","white","#f76a8c"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
               cluster_cols=FALSE,show_colnames = T,scale = "row", width = 1.8, height=5.5, legend = FALSE, annotation_legend =FALSE,annotation_row = NA,
               # annotation_row=absolute.promoter.activity.mean.tissue.cluster1[, "category", drop=FALSE], 
               # annotation_colors = list(category=color.selected),
               filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.heatmap.",cluster1.name,".pdf"))
               )
      
      ##Plot corresponding gene expression separately in a violin plot
      gene.expression.star.mean.tissue = gene.expression.star.mean[cluster1.genename,grep(name, colnames(gene.expression.star.mean))]
      colnames(gene.expression.star.mean.tissue) = sapply(strsplit(colnames(gene.expression.star.mean.tissue),"_"), "[", 2)
      gene.expression.star.mean.tissue.plot = melt(gene.expression.star.mean.tissue)
      gene.expression.star.mean.tissue.plot$Var2 = factor(gene.expression.star.mean.tissue.plot$Var2, levels = stage.order)
      line.plot = significant.isoform.tissue.expression.plot[significant.isoform.tissue.expression.plot$Type=="Gene",]
      
      p2 = ggplot(gene.expression.star.mean.tissue.plot, aes(x=Var2, y=value), group=1) +
        # geom_violin(alpha = 0.5, color="#69DADB", fill="#69DADB") +
        geom_boxplot(alpha = 0.5, outlier.color = "lightgrey",outlier.shape = NA, width = 0.3, color = "black") +  #
        # stat_compare_means(label = "p.signif",method = "wilcox.test",ref.group = stage.order[1],na.rm = TRUE,size=rel(2), hide.ns=TRUE, label.y = 23) + #length(stage.order)
        geom_line(data=line.plot, aes(x=stage, y=mean, group = 1), color = "black") +
        xlab("") +
        ylab("") +
        scale_y_continuous(limits = c(0,25)) +
        # ylab("log2(Gene Counts)")+
        theme_classic.add + theme(legend.position="none", axis.text.x = element_blank())
      print(p2)
      ggsave(p2, file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.genevioplot.",cluster1.name,".pdf")), width=2.2, height=2.3) #Fig. 2a-h, S2-7 upper pannel
      
      ##Plot corresponding gene heatmap
      pheatmap(gene.expression.star.mean.tissue, 
               color = colorRampPalette(c("#77d8d8","white","#f76a8c"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
               cluster_cols=FALSE,show_colnames = T,scale = "row", width = 1.8, height=5.5, legend = FALSE, annotation_legend =FALSE,annotation_row = NA,
               # annotation_row=absolute.promoter.activity.mean.tissue.cluster1[, "category", drop=FALSE], 
               filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.geneheatmap.",cluster1.name,".pdf")) #Fig. 2a-h, S2-7 bottom pannel
               )
    }  
    
    ##Type 1
    ##co-up and co-down regulated
    build.cluster.plot(cluster1.name="Up_Up", cluster1 = expand.grid(up.cluster, up.cluster))
    build.cluster.plot(cluster1.name="Down_Down", cluster1 = expand.grid(down.cluster, down.cluster))
    
    ##Type 2
    ##up-down and down-up mixed regulated
    minor.name = colnames(significant.isoform.tissue.table$IsoTable)
    major.up.cluster = up.cluster
    minor.down.cluster = grep(as.character(down.cluster), grep(as.character(flat.cluster),minor.name,invert=TRUE,value=TRUE), value=TRUE)
    major.down.cluster = down.cluster
    minor.up.cluster = grep(as.character(up.cluster), grep(as.character(flat.cluster),minor.name,invert=TRUE,value=TRUE), value=TRUE)
    build.cluster.plot(cluster1.name="Up_Down", cluster1 = expand.grid(major.up.cluster, minor.down.cluster))
    build.cluster.plot(cluster1.name="Down_Up", cluster1 = expand.grid(major.down.cluster, minor.up.cluster))
    
    ##Type 3
    minor.name = colnames(significant.isoform.tissue.table$IsoTable)
    major.flat.down.cluster = grep(as.character(down.cluster), grep(as.character(up.cluster),minor.name,invert=TRUE,value=TRUE), value=TRUE)
    major.flat.up.cluster = grep(as.character(up.cluster), grep(as.character(down.cluster),minor.name,invert=TRUE,value=TRUE), value=TRUE)
    ##flat-up and flat-down (major-minor)
    build.cluster.plot(cluster1.name="Flat_Up", cluster1 = expand.grid(flat.cluster, major.flat.up.cluster))
    build.cluster.plot(cluster1.name="Flat_Down", cluster1 = expand.grid(flat.cluster, major.flat.down.cluster))
    ##up-flat and down-flat (major-minor)
    build.cluster.plot(cluster1.name="Up_Flat", cluster1 = expand.grid(up.cluster, flat.cluster))
    build.cluster.plot(cluster1.name="Down_Flat", cluster1 = expand.grid(down.cluster, flat.cluster))
    
    ######catergories of major, minor promoter pairs into one dataframe
    category=c("Up_Up","Down_Down","Up_Down","Down_Up","Flat_Up","Flat_Down","Up_Flat","Down_Flat")
    significant.isoform.tissue.table.category = rbind()
    for (cluster.name in category){
      major.minor.table = read.delim(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.",cluster.name,".txt")), row.names=1)
      major.minor.table$category = cluster.name
      significant.isoform.tissue.table.category = rbind(significant.isoform.tissue.table.category, major.minor.table)
    }
    save(significant.isoform.tissue.table.category, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.category.RData")))
    
    
    rm(list=setdiff(ls(),basic.objects))
    gc()
  }
  
}  

#########################################################
######catergories number statistics of developmental dynamic promoters (sankey plot like circos plot) Fig. 2i-k
#########################################################
  
for (name in tissue.order){#"Brain" c("Testis")
  message(name)
  if(TRUE){
    # name="Brain"
    
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("isoform.model.tissue.RData")))
    # load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see0.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster0.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.category.RData")))
    
    ###developmental dynamic isoforms category classification depending on isoform number per gene
    #barplot
    significant.isoform.tissue.plot = as.data.frame(table(significant.isoform.tissue.see0$NumIso.by.gene)) # each gene may contain 1 to several DETs, this is the statistics of each gene"s DET number
    p1=ggplot(data = significant.isoform.tissue.plot[,], aes(x=Var1, y=Freq, width = 0.8)) +
      geom_bar(position="stack",stat="identity", fill = "#de7119",color = "#de7119", alpha = 0.6) + 
      labs(y = "Numbers of genes", x= "Promoter Number Per Gene") + 
      coord_flip() +
      # scale_x_discrete(c("Up","Down","Mixed")) +
      # scale_fill_manual(values =c("#ff1e56", "#ffac41", "#323232")) +
      # scale_color_manual(values =c("#ff1e56", "#ffac41", "#323232")) +
      ggtitle(name) +
      theme_classic.add + 
      theme(axis.title.x=element_blank(), #axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    p1
    ggsave(p1,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("seegenes.det.number.barplot.pdf")),width=2.5,height = 2.5)
   
    ##overlap between up and down clusters from significant.isoform.tissue.see between major and minor isoforms
    significant.isoform.tissue.table.plot = as.data.frame(table(significant.isoform.tissue.table.category$category))
    significant.isoform.tissue.table.plot = significant.isoform.tissue.table.plot[order(significant.isoform.tissue.table.plot$Var1),]
    #sankey plot
    sig.plot = significant.isoform.tissue.table.plot
    plot.color = c("#f76a8c", "#77d8d8", "#888888", "#f76a8c", "#77d8d8", "#888888") ##correspond to plot.node label
    # plot.color = c("#FB7AFC", "#77d8d8", "#888888", "#FB7AFC", "#77d8d8", "#888888") ##correspond to plot.node label
    plot.node = list(
      label = c("Major Up", "Major Down","Major Flat", "Minor Up", "Minor Down","Minor Flat"), 
      #############0            1               2               3           4             5
      ##labels of sources and targets c(0,1,2,3,4,5) (bars), correspond to bar number sig.plot$Var1
      color = plot.color,pad = 15,thickness = 20,
      line = list(color = "white", width = 0.5)
    )
    plot.link = list(
      source = c(1,1,1,2,2,0,0,0), ##According to sig.plot$Var1
      target = c(4,5,3,4,3,4,5,3),
      value = sig.plot$Freq ,
      color = c("rgba(119, 216, 216,0.3)", "rgba(119, 216, 216,0.3)", "rgba(119, 216, 216,0.3)", "rgba(136, 136, 136,0.3)", "rgba(136, 136, 136,0.3)", "rgba(247, 106, 140, 0.3)","rgba(247, 106, 140, 0.3)", "rgba(247, 106, 140,0.3)") #correspond to source/target number
    )
    
    fig <- plot_ly( type = "sankey", alpha = 0.3, orientation = "h", node = plot.node, link = plot.link)
    
    fig
    export(fig, file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("seegenes.cluster.5types.sankyplot.png")))
    
    #heatmap
    cluster.label = setNames(c("up", "down", "flat"), c(up.cluster,down.cluster,flat.cluster))
    significant.isoform.tissue.table.plot$Cluster.Mayor.name = gsub("_.*","",significant.isoform.tissue.table.plot$Var1)
    significant.isoform.tissue.table.plot$Cluster.minor.name = gsub(".*_","",significant.isoform.tissue.table.plot$Var1)
    p3=ggplot(data = significant.isoform.tissue.table.plot, aes(Cluster.Mayor.name, Cluster.minor.name, fill = log2(Freq)))+
      geom_tile(color = "white")+ 
      # geom_vline(xintercept = 3.5, linetype="dotted")+
      # geom_hline(yintercept = 3.5, linetype="dotted")+
      scale_fill_gradient2(low = "#00bcd4", high = "#FB7AFC",
                           na.value = "white", space = "Lab",
                           midpoint = median(log2(significant.isoform.tissue.table.plot$Freq)),
                           name="log2 (Overlap NO.)") +
      # scale_x_discrete(limits=as.character(c(1,4,5,2,3,6))) + 
      # scale_y_discrete(limits=as.character(c(1,4,5,2,3,6))) + 
      theme_minimal()+ theme_classic.add+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.3, size = 5, hjust = 0.3))+
      theme(axis.text.y = element_text(size = 5, hjust = 0.3))+
      labs(x=" ", y = " ")+
      coord_fixed() 
    p3
    ggsave(p3,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("seegenes.cluster.5types.heatmap.pdf")), width=3, height=3)
    
    rm(list=setdiff(ls(),basic.objects))
    gc()
  }
}


##############################################################
######promoter usage plot for co-up and co-down isoforms 
##############################################################
for (name in tissue.order){#"Brain" c("Testis")
  message(name)
  
  if(TRUE){
    # name = "Testis" #"Brain"
    
    ##load relative.promoter.activity info
    load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"relative.promoter.activity.star.mean.RData"))
    relative.promoter.activity.mean.tissue = relative.promoter.activity.star.mean[,grep(name, colnames(relative.promoter.activity.star.mean))]
    colnames(relative.promoter.activity.mean.tissue) = sapply(strsplit(colnames(relative.promoter.activity.mean.tissue), "_"), "[", 2)
    
    ##load up down cluster table category info
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.category.RData")))

    build.boxplot = function(data, promoter.type = "Major", file.name="up"){
      p<-ggplot(data[!is.na(data$value),], #subset(!is.na(relative.promoter.activity.mean.tissue.up.major)), 
                aes(x=Var2, y=value, color=color,fill=color)) +
        geom_boxplot(alpha=0.7,outlier.shape = NA, notch =FALSE, width=0.5) +
        scale_fill_manual(values =c("darkgrey", "#FF0075")) + #, "#323232" #ff1e56 #FE91CA #1DB9C3
        scale_color_manual(values =rep("black",2)) +
        scale_y_continuous(limits = quantile(data$value, c(0.2, 0.8),na.rm=TRUE)) + #c(0,1)
        labs(x="", y = "Relative Promoter Activity") +
        ggtitle(paste(name,promoter.type,"Promoter")) + theme_classic.add + 
        theme(plot.title = element_text(hjust = 0.5),legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
      if(file.name == "up"){
        p = p + stat_compare_means(label = "p.signif",method = "wilcox.test",hide.ns = TRUE, ref.group = as.character(data$Var2[1]),na.rm = TRUE,size=rel(2))
      } else if (file.name == "down"){
        p = p + stat_compare_means(label = "p.signif",method = "wilcox.test",hide.ns = TRUE, ref.group = as.character(data$Var2[length(as.character(data$Var2))-1]),na.rm = TRUE,size=rel(2))
      }
      print(p)
      ggsave(p, file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("promoterusage.boxplot.",file.name, promoter.type,".pdf")), width=3, height=3) #Fig. 5
    }
    
    ##up regulated cluster major and minor relative promoter activity along all stages
    up.cluster.table = significant.isoform.tissue.table.category[grep("Up", significant.isoform.tissue.table.category$category),]
    if(length(grep("Down", up.cluster.table$category))!=0){up.cluster.table = up.cluster.table[-grep("Down", up.cluster.table$category),]}
    relative.promoter.activity.mean.tissue.up.major = melt(relative.promoter.activity.mean.tissue[as.character(up.cluster.table$Major),])
    relative.promoter.activity.mean.tissue.up.major$color = gsub("[0-9].*","",relative.promoter.activity.mean.tissue.up.major$Var2)
    build.boxplot(data = relative.promoter.activity.mean.tissue.up.major, promoter.type = "Major", file.name="up")
    
    relative.promoter.activity.mean.tissue.up.minor = melt(relative.promoter.activity.mean.tissue[unlist(strsplit(as.character(up.cluster.table$Minor), "_")),])
    relative.promoter.activity.mean.tissue.up.minor$color = gsub("[0-9].*","",relative.promoter.activity.mean.tissue.up.minor$Var2)
    build.boxplot(data = relative.promoter.activity.mean.tissue.up.minor, promoter.type = "Minor", file.name="up")
    
    ##down regulated cluster major and minor relative promoter activity along all stages
    down.cluster.table = significant.isoform.tissue.table.category[grep("Down_Down", significant.isoform.tissue.table.category$category),]
    if(length(grep("Up", up.cluster.table$category))!=0){up.cluster.table = up.cluster.table[-grep("Up", up.cluster.table$category),]}
    relative.promoter.activity.mean.tissue.down.major = melt(relative.promoter.activity.mean.tissue[as.character(down.cluster.table$Major),])
    relative.promoter.activity.mean.tissue.down.major$color = gsub("[0-9].*","",relative.promoter.activity.mean.tissue.down.major$Var2)
    build.boxplot(data = relative.promoter.activity.mean.tissue.down.major, promoter.type = "Major", file.name="down")
    
    relative.promoter.activity.mean.tissue.down.minor = melt(relative.promoter.activity.mean.tissue[unlist(strsplit(as.character(down.cluster.table$Minor), "_")),])
    relative.promoter.activity.mean.tissue.down.minor$color = gsub("[0-9].*","",relative.promoter.activity.mean.tissue.down.minor$Var2)
    build.boxplot(data = relative.promoter.activity.mean.tissue.down.minor, promoter.type = "Minor", file.name="down")
    
    rm(list=setdiff(ls(),basic.objects))
    gc()
  }
}  
  

##############################################################
######promoter activity comparison between DDPs and unchanged promoters (Revision)
##############################################################
for (name in tissue.order){#"Brain" c("Testis")
  message(name)
  if(name=="Ovary"){gender="Female"} else {gender="Male"}
  
  if(TRUE){
    ##load the significant.isoform.tissue and expression info
    load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityCategory.star.RData"))
    load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolute.promoter.activity.star.mean.RData"))  
    absolute.promoter.activity.star.mean = absolute.promoter.activity.star.mean[,grep(gender, colnames(absolute.promoter.activity.star.mean))]
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.category.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster0.RData")))
  
    ##extract up,down minor promoters to unchanged promoters in the changed genes
    minor.promoters.updown = rownames(significant.isoform.tissue.see.cluster)
    minor.promoters.flat = rownames(significant.isoform.tissue.see.cluster0)
    
    # # Categorize promoters as changing or unchanging
    # minor_promoters <- absolutePromoterActivityCategory.star[absolutePromoterActivityCategory.star$PromoterCategory == "Minor", "promoterId"]
    # table(absolutePromoterActivityCategory.star[absolutePromoterActivityCategory.star$promoterId %in% minor.promoters.updown, "PromoterCategory"])
    # changing_promoters <- intersect(minor_promoters, minor.promoters.updown)
    # unchanging_promoters <- setdiff(minor_promoters, minor.promoters.updown)
    
    # Calculate the mean or max value for each promoter row
    # Use either rowMeans or apply for max
    promoter_values <- data.frame(
      promoterId = rownames(absolute.promoter.activity.star.mean),
      MaxValue = apply(absolute.promoter.activity.star.mean[,grep(name, colnames(absolute.promoter.activity.star.mean))], 1, max),
      MinValue = apply(absolute.promoter.activity.star.mean[,grep(name, colnames(absolute.promoter.activity.star.mean))], 1, min),
      MedianValue = apply(absolute.promoter.activity.star.mean[,grep(name, colnames(absolute.promoter.activity.star.mean))], 1, median),
      MeanValue = rowMeans(absolute.promoter.activity.star.mean[,grep(name, colnames(absolute.promoter.activity.star.mean))])
    )
    
    # Filter for minor promoters and merge with categories
    minor_promoter_values <- promoter_values[promoter_values$promoterId %in% c(minor.promoters.flat, minor.promoters.updown), ]
    # minor_promoter_values = minor_promoter_values[minor_promoter_values$meanValue != 0, ]
    minor_promoter_values$Category <- ifelse(minor_promoter_values$promoterId %in% minor.promoters.updown, "DDPs", "Flat_promoters")
    
    # Define colors for the categories
    colors <- c("DDPs" = "#de7119", "Flat_promoters" = "darkgrey")
    value_types <- c("MaxValue", "MeanValue","MinValue", "MedianValue")
    
    # Loop through each value type and generate/save the plot
    for (value_type in value_types) {
      
      # Dynamically set the y-axis label and title
      y_label <- paste(gsub("Value"," Value", value_type) ,"of Promoter Activity")
      plot_title <- ""
      
      # Create the boxplot
      plot <- ggplot(minor_promoter_values, aes(x = Category, y = .data[[value_type]], fill = Category, color = Category)) + 
        geom_boxplot(notch = TRUE, outlier.shape = NA, alpha=0.6) +  # Add notch and hide outliers
        scale_fill_manual(values = colors) +  # Apply custom colors
        scale_color_manual(values = colors) +  # Apply custom colors
        labs(
          title = plot_title,
          x = "",
          y = y_label
        ) +
        stat_compare_means(
          method = "wilcox.test",
          label = "p.format",
          na.rm = TRUE,size=rel(2.8), hide.ns=TRUE,
          ref.group = "DDPs"  # Compare the two groups directly
        ) +
        theme_classic.add  # Hide legend for simplicity
      
      # Save the plot
      file_name <- file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("promoteractivitycomparison.boxplot.",gsub("^M","m",value_type),".pdf")) #Fig. S8
      ggsave(file_name, plot=plot, width = 2.5, height = 2.5)    
    }
  }
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}


################################################
##Supplementary table S1 in the paper (Revision)
################################################
if(TRUE){
  promoter_annotation=rtracklayer::import(file.path(prewd,"ProcessedData","AnnotationObjects","promoter.bed"))
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityCategory.star.RData"))
  
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolute.promoter.activity.star.mean.RData"))  
  
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.category.RData")))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster.RData")))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster0.RData")))
  
  # Subset promoter_annotation to include only matching promoterIds
  promoter_annotation$promoterId <- gsub("_","",paste0("prmtr.", promoter_annotation$name))
  promoter_annotation$score <- NULL
  promoter_annotation$name <- NULL
  promoter_annotation_df <- as.data.frame(promoter_annotation)

  # Filter promoter_annotation_df to include only matching promoterIds
  filtered_promoter_annotation <- promoter_annotation_df %>%
    filter(promoterId %in% absolutePromoterActivityCategory.star$promoterId)
  
  # promoter activity in all tissues and stages
  absolute.promoter.activity.star.mean = as.data.frame(absolute.promoter.activity.star.mean)
  selected_columns <- c(
    grep("Male", colnames(absolute.promoter.activity.star.mean), value = TRUE),
    grep("Ovary", colnames(absolute.promoter.activity.star.mean), value = TRUE)
  )
  absolute.promoter.activity.star.mean = absolute.promoter.activity.star.mean[,selected_columns]
  colnames(absolute.promoter.activity.star.mean) = gsub("_Male","",gsub("_Female","",colnames(absolute.promoter.activity.star.mean)))
  absolute.promoter.activity.star.mean$promoterId <- rownames(absolute.promoter.activity.star.mean)
  
  # Merge the promoter granges with the category
  absolutePromoterActivityCategory.star_with_coords <- filtered_promoter_annotation %>% select(promoterId, seqnames, start, end, strand) %>% left_join( absolutePromoterActivityCategory.star, by = "promoterId")
  
  # Merge the promoter granges with the category, then merge with the promoter activity in all tissues, stages
  absolutePromoterActivityTable <- left_join(absolutePromoterActivityCategory.star_with_coords, absolute.promoter.activity.star.mean, by = "promoterId")  

  head(absolutePromoterActivityTable)
  save(absolutePromoterActivityTable, file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityTable.RData"))  
  write.table(absolutePromoterActivityTable, file=file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityTable.txt"), quote = FALSE, row.names = FALSE)  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
  
}


################################################
##Plot for each gene and isoform changes
################################################
  
for (name in tissue.order){#"Brain" c("Testis")
  message(name)
  
  if(TRUE){
    # name = "Brain"
    
    ##Plot DE isoforms found without control
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.majoriso.change.RData")))
    build.iso.plot = function(genename){
      pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.",genename,".pdf")))
      IsoPlot(significant.isoform.tissue,genename,only.sig.iso=FALSE,cex.main=2,cex.legend=1)
      dev.off()
    }
    registerDoParallel(cores = 10)
    foreach(n=significant.isoform.tissue.majoriso.change$L) %dopar% build.iso.plot(n)
    
    
    rm(list=setdiff(ls(),basic.objects))
    gc()
  }
  
  if(TRUE){
    ##intersect single tissue DE isoforms with significantly compared to other tissue
    ##here, the aim is to find DE isoforms only existed in this tissue while not others, or the changing pattern is different from the else
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.all.RData")))
    significant.isoform.merge.list = c(Testis=list(significant.isoform.tissue$DET), lapply(significant.isoform.all$summary[-1], as.character))
    significant.isoform.merge = Reduce(intersect, significant.isoform.merge.list )
    significant.isoform.merge.expression = significant.isoform.all$sig.genes[[2]]$sig.profiles[significant.isoform.merge,]
    
    ##tissue specific isoforms list in each tissue at a stage
    load(file.path(prewd,"ProcessedData","TissueObjects",method,"deseq2",paste0("promoter.tissue.specific.list.", "s9wpb_Male", ".RData") ))
    tissue.specific.list = promoter.tissue.specific.list[[name]]
    
    
  }
}

################################################
######profileplot and heatmap for 5 types isoforms corresponding genes
################################################
# if(TRUE){
#   # name = "Brain"
#   
#   ##Extract genes with up or down clusters in both major and minor
#   # sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects","human.sample.sheet.txt"),header=TRUE, stringsAsFactors = FALSE)
#   # sample.info = sample.info[sample.info$tissue=="Brain",]
#   # load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
#   cluster1 = read.delim(file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.up.txt")), stringsAsFactors = FALSE)
#   cluster2 = read.delim(file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.down.txt")), stringsAsFactors = FALSE)
#   
#   # ##Plot Major and minor profile in one plot
#   # load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("isoform.model.tissue.RData")))
#   # significant.isoform.tissue.selected.gene = significant.isoform.tissue$get2$sig.genes$gene$sig.profiles
#   # significant.isoform.tissue.selected.cluster1 = significant.isoform.tissue.selected.gene[rownames(significant.isoform.tissue.selected.gene) %in% unique(cluster1$X),]
#   # significant.isoform.tissue.selected.cluster2 = significant.isoform.tissue.selected.gene[rownames(significant.isoform.tissue.selected.gene) %in% unique(cluster2$X),]
#   # significant.isoform.tissue.selected = significant.isoform.tissue$get2
#   # 
#   # pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.cluster.up.gene.pdf")))
#   # PlotGroups(data = significant.isoform.tissue.selected.cluster1,
#   #               x.labels = unique(sample.info$stage),
#   #               edesign = significant.isoform.tissue.selected$edesign,
#   #               groups = significant.isoform.tissue.selected$edesign[, c(3:ncol(significant.isoform.tissue.selected$edesign)),drop=FALSE],
#   #               summary.mode = "mean",
#   #               groups.vector = significant.isoform.tissue.selected$groups.vector)
#   # dev.off()
#   # pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.cluster.down.gene.pdf")))
#   # PlotTwoGroups(data = significant.isoform.tissue.selected.cluster2,
#   #               data2 = significant.isoform.tissue.selected.minor.cluster2,
#   #               x.labels = unique(sample.info$stage),
#   #               edesign = significant.isoform.tissue.selected$edesign,
#   #               groups = significant.isoform.tissue.selected$edesign[, c(3:ncol(significant.isoform.tissue.selected$edesign)),drop=FALSE],
#   #               summary.mode = "mean",
#   #               groups.vector = significant.isoform.tissue.selected$groups.vector)
#   # dev.off()
#   # 
#   ##Plot Major and minor heatmap
#   load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
#   geneExpression.star = get(load(file.path(prewd,"ProcessedData","PromoterObjects","gene.expression.mean.RData")))
#   absolute.promoter.activity.mean.tissue = geneExpression.star[,grep(name,colnames(geneExpression.star))]
#   absolute.promoter.activity.mean.tissue.cluster1 = as.data.frame(absolute.promoter.activity.mean.tissue[unique(cluster1$X),]) #
#   absolute.promoter.activity.mean.tissue.cluster2 = as.data.frame(absolute.promoter.activity.mean.tissue[unique(cluster2$X),]) #
#   
#   pheatmap(absolute.promoter.activity.mean.tissue.cluster1, 
#            color = colorRampPalette(c("#00bcd4","white","#FB7AFC"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
#            cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=10, legend = TRUE,
#            filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.heatmap.up.gene.pdf")))
#   #annotation_colors = as.list(absolute.promoter.activity.mean.tissue.cluster1[,"color"]),
#   pheatmap(absolute.promoter.activity.mean.tissue.cluster2, 
#            color = colorRampPalette(c("#00bcd4","white","#FB7AFC"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
#            cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=10, legend = TRUE,
#            filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.heatmap.down.gene.pdf")))
#   
#   
#   ######profileplot and heatmap for different trend isoforms 
#   ################################################
#   cluster1 = read.delim(file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.mixedup.txt")), stringsAsFactors = FALSE)
#   cluster2 = read.delim(file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.mixeddown.txt")), stringsAsFactors = FALSE)
#   
#   # ##Major and minor promter name belong to extracted genes
#   # major.cluster1 = as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,"Major"])
#   # major.cluster2 = as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,"Major"])
#   # build.strsplit = function(x){return(unique(unlist(strsplit(x,"_"))))}
#   # minor.cluster1 = build.strsplit(as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,"Minor"]))
#   # minor.cluster2 = build.strsplit(as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,"Minor"]))
#   # write.table(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,], file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.mixedup.txt")), sep="\t", col.names = NA, row.names = TRUE, quote = FALSE)
#   # write.table(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,], file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.mixeddown.txt")), sep="\t", col.names = NA, row.names = TRUE, quote = FALSE)
#   # 
#   # ##Plot Major and minor profile in one plot
#   # load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
#   # significant.isoform.tissue.selected = significant.isoform.tissue$get2$sig.genes
#   # significant.isoform.tissue.selected.minor.cluster1 = significant.isoform.tissue.selected[rownames(significant.isoform.tissue.selected) %in% minor.cluster1,]
#   # significant.isoform.tissue.selected.minor.cluster2 = significant.isoform.tissue.selected[rownames(significant.isoform.tissue.selected) %in% minor.cluster2,]
#   # significant.isoform.tissue.selected.major.cluster1 = significant.isoform.tissue.selected[rownames(significant.isoform.tissue.selected) %in% major.cluster1,]
#   # significant.isoform.tissue.selected.major.cluster2 = significant.isoform.tissue.selected[rownames(significant.isoform.tissue.selected) %in% major.cluster2,]
#   # 
#   # pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.cluster.mixedup.gene.pdf")))
#   # PlotTwoGroups(data = significant.isoform.tissue.selected.major.cluster1, 
#   #               data2 = significant.isoform.tissue.selected.minor.cluster1, 
#   #               x.labels = unique(sample.info$stage),
#   #               edesign = significant.isoform.tissue.selected$edesign, 
#   #               groups = significant.isoform.tissue.selected$edesign[, c(3:ncol(significant.isoform.tissue.selected$edesign)),drop=FALSE], 
#   #               summary.mode = "mean",
#   #               groups.vector = significant.isoform.tissue.selected$groups.vector)
#   # dev.off()
#   # pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.cluster.mixeddown.gene.pdf")))
#   # PlotTwoGroups(data = significant.isoform.tissue.selected.major.cluster2, 
#   #               data2 = significant.isoform.tissue.selected.minor.cluster2,
#   #               x.labels = unique(sample.info$stage),
#   #               edesign = significant.isoform.tissue.selected$edesign, 
#   #               groups = significant.isoform.tissue.selected$edesign[, c(3:ncol(significant.isoform.tissue.selected$edesign)),drop=FALSE], 
#   #               summary.mode = "mean",
#   #               groups.vector = significant.isoform.tissue.selected$groups.vector)
#   # dev.off()
#   
#   ##Plot Major and minor heatmap
#   geneExpression.star = get(load(file.path(prewd,"ProcessedData","PromoterObjects","gene.expression.mean.RData")))
#   absolute.promoter.activity.mean.tissue = geneExpression.star[,grep(name,colnames(geneExpression.star))]
#   absolute.promoter.activity.mean.tissue.cluster1 = as.data.frame(absolute.promoter.activity.mean.tissue[unique(cluster1$X),]) #
#   absolute.promoter.activity.mean.tissue.cluster2 = as.data.frame(absolute.promoter.activity.mean.tissue[unique(cluster2$X),]) #
#   
#   pheatmap(absolute.promoter.activity.mean.tissue.cluster1, 
#            color = colorRampPalette(c("#00bcd4","white","#FB7AFC"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
#            cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=10, legend = TRUE,
#            filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.heatmap.mixedup.gene.pdf")))
#   #annotation_colors = as.list(absolute.promoter.activity.mean.tissue.cluster1[,"color"]),
#   pheatmap(absolute.promoter.activity.mean.tissue.cluster2, 
#            color = colorRampPalette(c("#00bcd4","white","#FB7AFC"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
#            cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=10, legend = TRUE,
#            filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.heatmap.mixeddown.gene.pdf")))
#   
#   rm(list=setdiff(ls(),basic.objects))
#   gc()
# }




# ################################################
# ###With control, Time series analysis for different developmental stages in each tissue
# ################################################
# if(TRUE){
#   dir.create(file.path(prewd,"ProcessedData","ControlPromoterObjects"))
#   dir.create(file.path(prewd,"Output","ControlStageSpecific"))
#   
#   ###different promoters contributed difference along the development
#   ##load information
#   sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects","human.sample.sheet.txt"),header=TRUE, stringsAsFactors = FALSE)
#   absolutePromoterActivity.star = get(load(file.path(prewd,"ProcessedData","PromoterObjects","absolutePromoterActivity.star.cpm.RData")))
#   # load(file.path(prewd,"ProcessedData","PromoterObjects","absolutePromoterActivity.star.RData"))
#   absolutePromoterActivity.star[is.na(absolutePromoterActivity.star)] = 0
#   
#   geneExpression.star = get(load(file.path(prewd,"ProcessedData","PromoterObjects","geneExpression.star.cpm.RData")))
#   # load(file.path(prewd,"ProcessedData","PromoterObjects","geneExpression.star.RData"))
#   load(file.path(prewd,"ProcessedData","AnnotationObjects","promoterAnnotationData.RData"))
#   geneExpression.star = merge(geneExpression.star,unique(promoterIdMapping(promoterAnnotationData)[,c("promoterId","geneId")]), by.x="geneId", by.y="geneId")
#   rownames(geneExpression.star) = geneExpression.star$promoterId
#   geneExpression.star = geneExpression.star[,-ncol(geneExpression.star)]
#   
#   ##With Control (Gene expression combined all isoforms)
#   ####################################
#   ##Counts information  in each tissue
#   name = "Brain"
#   sample.info.tissue = sample.info[sample.info$tissue %in% name,] #,"Cerebellum"
#   dir.create(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control")))
#   dir.create(file.path(prewd,"Output","StageSpecific",method,norm.method,paste0(name,"Control")))
#   
#   ##normalized counts is used for promoter in each tissue (Treatment)
#   absolutePromoterActivity.tissue = absolutePromoterActivity.star[,sample.info.tissue$samplename]
#   absolutePromoterActivity.tissue = 2^(absolutePromoterActivity.tissue)
#   absolutePromoterActivity.tissue[absolutePromoterActivity.tissue==1] = 0
#   absolutePromoterActivity.tissue = absolutePromoterActivity.tissue[rowSums(absolutePromoterActivity.tissue)>=10 & 
#                                                                       apply(absolutePromoterActivity.tissue,1,function(x){max(x)>1}),]
#   # absolutePromoterActivity.tissue = cpm(absolutePromoterActivity.tissue, lib.size = NULL)
#   # absolutePromoterActivity.tissue.control = do.call(cbind.data.frame, rep(absolutePromoterActivity.tissue[,2,drop=FALSE], ncol(absolutePromoterActivity.tissue)-1))
#   
#   ##normalized counts is used for correspond gene in each tissue (Control)
#   absolutePromoterActivity.tissue.control = geneExpression.star[rownames(absolutePromoterActivity.tissue),colnames(absolutePromoterActivity.tissue)]
#   colnames(absolutePromoterActivity.tissue.control) = paste("Control", sapply(strsplit(colnames(absolutePromoterActivity.tissue.control), "_"), function(x){return(paste(x[-1], collapse = "_"))}), sep="_")
#   absolutePromoterActivity.tissue.control = 2^(absolutePromoterActivity.tissue.control)
#   absolutePromoterActivity.tissue.control[absolutePromoterActivity.tissue.control==1] = 0
#   
#   absolutePromoterActivity.tissue.combine = cbind(absolutePromoterActivity.star[rownames(absolutePromoterActivity.tissue), "geneId", drop=FALSE],absolutePromoterActivity.tissue, absolutePromoterActivity.tissue.control)
#   
#   ####################################
#   ##sample info for the design
#   sample.info.tissue.control = sample.info.tissue
#   sample.info.tissue.control$tissue = "Control"
#   sample.info.tissue.control$condition = paste(sample.info.tissue.control$tissue, sample.info.tissue.control$stage, sep="_")
#   sample.info.tissue.control$samplename = paste("Control", sapply(strsplit(sample.info.tissue.control$samplename, "_"), function(x){return(paste(x[-1], collapse = "_"))}), sep="_")
#   sample.info.tissue.combine = rbind(sample.info.tissue, sample.info.tissue.control)
#   rownames(sample.info.tissue.combine) = sample.info.tissue.combine$samplename
#   
#   time.change = gsub("6month","1.5",gsub("0month", "0.76", gsub("year","",gsub("week", "", sample.info.tissue.combine$stage))))  #1:length(unique(sample.info.tissue.combine$stage))
#   time = 1:length(time.change)
#   time[grep("b",time.change)] = log2(as.numeric(gsub("b", "", time.change[grep("b",time.change)])) * 7)
#   time[grep("p",time.change)] = log2(as.numeric(gsub("p", "", time.change[grep("p",time.change)])) * 365)
# 
#   sample.info.tissue.combine$time = time #[sample.info.tissue.combine$stage]
#   replicate = 1:length(unique(sample.info.tissue.combine$condition)) #paste(sample.info.tissue.combine$stage,sample.info.tissue.combine$tissue,sep="_")
#   names(replicate) = unique(sample.info.tissue.combine$condition)
#   sample.info.tissue.combine$replicate = replicate[sample.info.tissue.combine$condition]  #unlist(sapply(split(sample.info.tissue.combine$stage,factor(replicate,levels = unique(replicate))),function(x)return(1:length(x))))
#   sample.info.tissue.combine$gene = as.numeric(sample.info.tissue.combine$tissue=="Control")
#   sample.info.tissue.combine$promoter = as.numeric(sample.info.tissue.combine$tissue!="Control")
#   sample.info.tissue.combine = sample.info.tissue.combine[,-1:-4]
#   
#   ####################################
#   ##regression fit
#   design.tissue <- make.design.matrix(sample.info.tissue.combine, degree = 3)
#   isoform.model.tissue <- IsoModel(data=absolutePromoterActivity.tissue.combine[,-1], gen=absolutePromoterActivity.tissue.combine[,1], design=design.tissue, Q=1, Qfit = 0.05, counts = TRUE) #, counts=FALSE, Q = 0.05, min.obs=15
#   save(isoform.model.tissue, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("isoform.model.tissue.RData")))
#   
#   significant.isoform.tissue <- getDS(Model = isoform.model.tissue, rsq = 0.3, vars = "all") #groups
#   length(significant.isoform.tissue$List0)
#   length(significant.isoform.tissue$DSG)
#   length(significant.isoform.tissue$DET)
#   save(significant.isoform.tissue, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("significant.isoform.tissue.RData")))
#   
#   ##clustering
#   pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,paste0(name,"Control"),paste0("seegenes.cluster.pdf")))
#   significant.isoform.tissue.see <- seeDS(significant.isoform.tissue, cluster.all=FALSE, k=9)
#   dev.off()
#   save(significant.isoform.tissue.see, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("significant.isoform.tissue.see.RData")))
#   
#   ##find isoforms belong to same gene but belong to different clusters
#   significant.isoform.tissue.table <- tableDS(significant.isoform.tissue.see)
#   significant.isoform.tissue.table$IsoTable
#   getDSPatterns(significant.isoform.tissue.table, 2, 5)
#   save(significant.isoform.tissue.table, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("significant.isoform.tissue.table.RData")))
#   
#   ##find differentially changed alternative promoters comparing gene changes to promoter changes
#   significant.isoform.tissue.majoriso.change = PodiumChange(get = significant.isoform.tissue, only.sig.iso=TRUE, comparison="any")
#   length(significant.isoform.tissue.majoriso.change$L)
#   save(significant.isoform.tissue.majoriso.change, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("significant.isoform.tissue.majoriso.change.RData")))
#   
#   
#   rm(list=setdiff(ls(),basic.objects))
#   gc()
#   
# }
# 
# 
# ################################################
# ##Plot for each gene and isoform changes
# ################################################
# 
# if(TRUE){
#   
#   ##Plot DE isoforms found without control
#   load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("significant.isoform.tissue.RData")))
#   load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,paste0("significant.isoform.tissue.majoriso.change.RData")))
#   build.iso.plot = function(genename){
#     pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,paste0("majoriso.change.",genename,".pdf")))
#     IsoPlot(significant.isoform.tissue,genename,only.sig.iso=FALSE,cex.main=2,cex.legend=1)
#     dev.off()
#   }
#   registerDoParallel(cores = 10)
#   foreach(n=significant.isoform.tissue.majoriso.change$L) %dopar% build.iso.plot(n)
#   
#   # ##Plot DE isoforms found with control
#   # load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("significant.isoform.tissue.RData")))
#   # load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,paste0(name,"Control"),paste0("significant.isoform.tissue.majoriso.change.RData")))
#   # build.iso.plot = function(genename){
#   #   pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,paste0(name,"Control"),paste0("majoriso.change.",genename,".pdf")))
#   #   IsoPlot(significant.isoform.tissue,genename,only.sig.iso=FALSE,cex.main=2,cex.legend=1)
#   #   dev.off()
#   # }
#   # registerDoParallel(cores = 10)
#   # foreach(n=significant.isoform.tissue.majoriso.change$L) %dopar% build.iso.plot(n)
#   
#   # ##Plot major and minor isoform clusters
#   # ##For major isoform, classify into two clusters (one is up and one is down)
#   # significant.promoter.major = as.character(significant.isoform.tissue.table$IsoMajorMinor[significant.isoform.tissue.majoriso.change$L,"Major"])
#   # significant.isoform.tissue.selected.major = significant.isoform.tissue$get2$sig.genes
#   # significant.isoform.tissue.selected.major$sig.profiles = significant.isoform.tissue.selected.major$sig.profiles[significant.promoter.major,]
#   # significant.isoform.tissue.selected.major$sig.pvalues = significant.isoform.tissue.selected.major$sig.pvalues[significant.promoter.major,]
#   # significant.isoform.tissue.selected.major$coefficients = significant.isoform.tissue.selected.major$coefficients[significant.promoter.major,]
#   # significant.isoform.tissue.selected.major$group.coeffs = significant.isoform.tissue.selected.major$group.coeffs[significant.promoter.major,]
#   # significant.isoform.tissue.selected.major$g = length(significant.promoter.major)
#   # pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,paste0("majoriso.change.cluster.major.pdf")))
#   # significant.isoform.tissue.selected.major.see = see.genes(significant.isoform.tissue.selected.major, item = "Isoforms", cluster.method="hclust" ,k = 2, newX11 = FALSE)
#   # dev.off()
#   # major.cluster1 = names(significant.isoform.tissue.selected.major.see$cut[significant.isoform.tissue.selected.major.see$cut == 1])
#   # major.cluster2 = names(significant.isoform.tissue.selected.major.see$cut[significant.isoform.tissue.selected.major.see$cut == 2])
#   
#   rm(list=setdiff(ls(),basic.objects))
#   gc()
# }


# ################################################
# ####Stage specific analysis
# ################################################
# ###development specific identification###
# if(TRUE){
#   sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects","human.sample.sheet.txt"),header=TRUE, stringsAsFactors = FALSE)
#   
#   load(file.path(prewd,"ProcessedData","PromoterObjects","absolute.promoter.activity.mean.RData"))
#   # absolute.promoter.activity.mean = absolute.promoter.activity.mean[,-grep("Forebrain",colnames(absolute.promoter.activity.mean))]
#   # absolute.promoter.activity.mean = absolute.promoter.activity.mean[,-grep("Hindbrain",colnames(absolute.promoter.activity.mean))]
#   
#   ###development specific score calculation in each developmental stage using tau and spm###
#   all.tissue = unique(sample.info$tissue)
#   absolute.promoter.activity.mean.stage.specific.tau = cbind()
#   for(sample in all.tissue){
#     absolute.promoter.activity.mean.stage.specific.tau = cbind(absolute.promoter.activity.mean.stage.specific.tau, apply(absolute.promoter.activity.mean[,grep(sample,colnames(absolute.promoter.activity.mean))], 1, fTau))
#   }
#   colnames(absolute.promoter.activity.mean.stage.specific.tau) = all.tissue
#   save(absolute.promoter.activity.mean.stage.specific.tau,file=file.path(prewd,"ProcessedData","PromoterObjects","absolute.promoter.activity.mean.stage.specific.tau.RData"))
#   
#   absolute.promoter.activity.mean.stage.specific.spm = cbind()
#   for(sample in all.tissue){
#     absolute.promoter.activity.mean.stage.specific.spm = cbind(absolute.promoter.activity.mean.stage.specific.spm, apply(absolute.promoter.activity.mean[,grep(sample,colnames(absolute.promoter.activity.mean))], 1, fSpm))
#   }
#   colnames(absolute.promoter.activity.mean.stage.specific.spm) = all.tissue
#   save(absolute.promoter.activity.mean.stage.specific.spm,file=file.path(prewd,"ProcessedData","PromoterObjects","absolute.promoter.activity.mean.stage.specific.spm.RData"))
#   
#   ##In each development-stage, compare with other development-stage,  mean absolute promoter activity development-stage/other>= 2; mean gene expression development-stage/other<= 1.5
#   load(file.path(prewd,"ProcessedData","PromoterObjects","gene.expression.mean.RData"))
#   gene.expression.mean = gene.expression.mean[,colnames(absolute.promoter.activity.mean)]
#   load(file.path(prewd,"ProcessedData","PromoterObjects","relative.promoter.activity.mean.RData"))
#   relative.promoter.activity.mean[is.na(relative.promoter.activity.mean)] = 0
#   relative.promoter.activity.mean = relative.promoter.activity.mean[,colnames(absolute.promoter.activity.mean)]
#   
#   load(file.path(prewd,"ProcessedData","AnnotationObjects","promoterAnnotationData.RData"))
#   load(file.path(prewd,"ProcessedData","PromoterObjects","absolute.promoter.activity.mean.stage.specific.tau.RData"))
#   load(file.path(prewd,"ProcessedData","PromoterObjects","absolute.promoter.activity.mean.stage.specific.spm.RData"))
#   
#   promoter.stage.specific.list=list()
#   for(sample in all.tissue){
#     absolute.promoter.activity.mean.selected = absolute.promoter.activity.mean[,grep(sample,colnames(absolute.promoter.activity.mean))]
#     relative.promoter.activity.mean.selected = relative.promoter.activity.mean[,grep(sample,colnames(relative.promoter.activity.mean))]
#     gene.expression.mean.selected = gene.expression.mean[,grep(sample,colnames(gene.expression.mean))]
#     
#     all.stage = unique(sapply((strsplit(colnames(absolute.promoter.activity.mean.selected),"_")),function(x){return(x[2])}))
#     for(stage in all.stage){
#       promoter.activity.mean.all = data.frame(absolute.stage = rowMeans(absolute.promoter.activity.mean.selected[,grep(stage,colnames(absolute.promoter.activity.mean.selected)),drop=FALSE]),
#                                               absolute.other = rowMeans(absolute.promoter.activity.mean.selected[,-grep(stage,colnames(absolute.promoter.activity.mean.selected)),drop=FALSE]),
#                                               relative.stage = rowMeans(relative.promoter.activity.mean.selected[,grep(stage,colnames(relative.promoter.activity.mean.selected)),drop=FALSE]),
#                                               relative.other = rowMeans(relative.promoter.activity.mean.selected[,-grep(stage,colnames(relative.promoter.activity.mean.selected)),drop=FALSE]),
#                                               tau = absolute.promoter.activity.mean.stage.specific.tau[rownames(absolute.promoter.activity.mean.selected),sample],
#                                               spm = absolute.promoter.activity.mean.stage.specific.spm[rownames(absolute.promoter.activity.mean.selected),sample],
#                                               stringsAsFactors = FALSE)
#       
#       gene.expression.mean.all = data.frame(expression.stage = rowMeans(gene.expression.mean.selected[,grep(stage,colnames(gene.expression.mean.selected)),drop=FALSE]),
#                                             expression.other = rowMeans(gene.expression.mean.selected[,-grep(stage,colnames(gene.expression.mean.selected)),drop=FALSE]),
#                                             stringsAsFactors = FALSE)
#       
#       promoter.activity.mean.all.merge = merge(promoter.activity.mean.all,unique(promoterAnnotationData@promoterIdMapping[,c("promoterId","geneId")]),by.x=0,by.y="promoterId")
#       promoter.activity.mean.all.merge = merge(promoter.activity.mean.all.merge,gene.expression.mean.all,by.x="geneId",by.y=0)
#       rownames(promoter.activity.mean.all.merge) = promoter.activity.mean.all.merge$Row.names
#       promoter.activity.mean.all.merge$absolute.divide = promoter.activity.mean.all.merge$absolute.stage/promoter.activity.mean.all.merge$absolute.other
#       promoter.activity.mean.all.merge$expression.divide = promoter.activity.mean.all.merge$expression.stage/promoter.activity.mean.all.merge$expression.other
#       promoter.activity.mean.all.merge[is.na(promoter.activity.mean.all.merge)] = 0
#       
#       promoter.activity.stage.specific = promoter.activity.mean.all.merge[promoter.activity.mean.all.merge$absolute.stage>=0.25 & 
#                                                                                   promoter.activity.mean.all.merge$absolute.other>=0.25 &
#                                                                                   promoter.activity.mean.all.merge$relative.stage>=0.25 & 
#                                                                                   promoter.activity.mean.all.merge$relative.other>=0.25 & 
#                                                                                   promoter.activity.mean.all.merge$expression.stage>=1 & 
#                                                                                   promoter.activity.mean.all.merge$expression.other>=1 & 
#                                                                                   promoter.activity.mean.all.merge$absolute.divide>=2 & 
#                                                                                   promoter.activity.mean.all.merge$expression.divide<=1.5 & 
#                                                                                   (promoter.activity.mean.all.merge$tau >=0.5 | promoter.activity.mean.all.merge$spm >=0.3),,drop=FALSE]
#       
#       promoter.stage.specific.list[[paste(sample,stage,sep="_")]] = promoter.activity.stage.specific
#     }
#     
#   }
#   save(promoter.stage.specific.list,file=file.path(prewd,"ProcessedData","PromoterObjects","promoter.stage.specific.list.RData"))
#   
#   
#   rm(list=setdiff(ls(),basic.objects))
#   gc()
# }

# #barplot
# p2=ggplot(data = significant.isoform.tissue.table.plot, aes(x=category, y=Freq, width = 0.8)) +
#   geom_bar(position="stack",stat="identity", fill = "#de7119",color = "#de7119", alpha = 0.6) + 
#   labs(y = "Numbers of genes") + 
#   scale_x_discrete(c("Up","Down","Mixed")) +
#   # scale_fill_manual(values =c("#ff1e56", "#ffac41", "#323232")) +
#   # scale_color_manual(values =c("#ff1e56", "#ffac41", "#323232")) +
#   ggtitle(name) +
#   theme_classic.add + 
#   theme(axis.title.x=element_blank(), #axis.text.x=element_blank(), axis.ticks.x=element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())
# p2
# ggsave(p2,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,paste0("seegenes.cluster.3types.barplot.pdf")),width=2.5,height = 4)


##one isoform is developmental dynamically changing while the others are not, extract and plot all of them
# promoter.gene.correspond = data.frame(promoterid = rownames(significant.isoform.tissue$Model$data), geneid = significant.isoform.tissue$Model$gen, stringsAsFactors = FALSE)
# significant.isoform.tissue.transcript.number.only1 = significant.isoform.tissue.see.cluster[significant.isoform.tissue.see.cluster$numiso==1,"geneid"]
# significant.isoform.tissue.transcript.number.only1.all.isoforms = promoter.gene.correspond[promoter.gene.correspond$geneid %in% significant.isoform.tissue.transcript.number.only1,]
# significant.isoform.tissue.transcript.number.only1.all.isoforms$cluster = 0
# significant.isoform.tissue.transcript.number.only1.all.isoforms$numiso = 1
# 
# significant.isoform.tissue.transcript.number.only1.up = significant.isoform.tissue.see.cluster[significant.isoform.tissue.see.cluster$numiso==1 & significant.isoform.tissue.see.cluster$cluster==up.cluster,"geneid"]
# significant.isoform.tissue.transcript.number.only1.down = significant.isoform.tissue.see.cluster[significant.isoform.tissue.see.cluster$numiso==1 & significant.isoform.tissue.see.cluster$cluster==down.cluster,"geneid"]

# significant.isoform.tissue.transcript.number.only1.all.isoforms.up = significant.isoform.tissue$Model$data[significant.isoform.tissue$Model$gen %in% significant.isoform.tissue.transcript.number.only1.up,]
# significant.isoform.tissue.transcript.number.only1.all.isoforms.down = significant.isoform.tissue$Model$data[significant.isoform.tissue$Model$gen %in% significant.isoform.tissue.transcript.number.only1.down,]
# significant.isoform.tissue.transcript.number.only1.all.isoforms = significant.isoform.tissue$Model$data[significant.isoform.tissue$Model$gen %in% significant.isoform.tissue.transcript.number.only1,]
# significant.isoform.tissue.transcript.number.only1.all.isoforms.left = significant.isoform.tissue.transcript.number.only1.all.isoforms[(!rownames(significant.isoform.tissue.transcript.number.only1.all.isoforms) %in% rownames(significant.isoform.tissue.see.cluster)),]

##for 1 isoform changed genes, how the other isoforms behave, so we clustered them, but not used anymore
# extract.name = rownames(significant.isoform.tissue.transcript.number.only1.all.isoforms)
# extract.model = isoform.model.tissue.all$Tfit #
# significant.isoform.tissue.inactive = significant.isoform.tissue$get2$sig.genes
# significant.isoform.tissue.inactive$sig.profiles = extract.model$dat[extract.name,]
# significant.isoform.tissue.inactive$coefficients = extract.model$coefficients[extract.name,]
# significant.isoform.tissue.inactive$group.coeffs = extract.model$group.coeffs[extract.name,]
# significant.isoform.tissue.inactive$g = dim(significant.isoform.tissue.inactive$sig.profiles)
# pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("seegenes.other.nosignificant.isoforms.pdf")))
# significant.isoform.tissue.inactive.see = see.genes(data = significant.isoform.tissue.inactive, item = "Isoforms", cluster.method = "hclust", k = 4, newX11 = FALSE)
# dev.off()
# save(significant.isoform.tissue.inactive.see, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.inactive.see.RData")))
# significant.isoform.tissue.inactive.table <- tableDS(significant.isoform.tissue.inactive.see)
# print(significant.isoform.tissue.inactive.table$IsoTable)

# if(TRUE){
# pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.lineplot.",cluster1.name,".pdf")))
# print(PlotTwoGroups(data = significant.isoform.tissue.expression.major.cluster1, 
#                     data2 = significant.isoform.tissue.expression.minor.cluster1, 
#                     x.labels = stage.info,
#                     edesign = significant.isoform.tissue.see0$Model$design$edesign, 
#                     groups = significant.isoform.tissue.see0$Model$design$edesign[, c(3:ncol(significant.isoform.tissue.see0$Model$design$edesign)),drop=FALSE], 
#                     summary.mode = "mean")) #, groups.vector = significant.isoform.tissue$groups.vector
# dev.off()
# 
# cluster2 = expand.grid(down.cluster, down.cluster) ##major down regulated
# cluster2.genename = build.extract.cluster.name(cluster2)
# # cluster2.genename = cluster2.genename[cluster2.genename %in% significant.isoform.tissue.majoriso.change$L]
# major.cluster2 = as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,"Major"])
# minor.cluster2 = build.strsplit(as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,"Minor"]))
# write.table(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,], file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.down.txt")), sep="\t", col.names = NA, row.names = TRUE, quote = FALSE)
# significant.isoform.tissue.selected.minor.cluster2 = significant.isoform.tissue.selected$sig.profiles[rownames(significant.isoform.tissue.selected$sig.profiles) %in% minor.cluster2,]
# significant.isoform.tissue.selected.major.cluster2 = significant.isoform.tissue.selected$sig.profiles[rownames(significant.isoform.tissue.selected$sig.profiles) %in% major.cluster2,]
# 
# pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.cluster.down.pdf")))
# PlotTwoGroups(data = significant.isoform.tissue.selected.major.cluster2, 
#               data2 = significant.isoform.tissue.selected.minor.cluster2,
#               x.labels = unique(sample.info$stage),
#               edesign = significant.isoform.tissue.selected$edesign, 
#               groups = significant.isoform.tissue.selected$edesign[, c(3:ncol(significant.isoform.tissue.selected$edesign)),drop=FALSE], 
#               summary.mode = "mean",
#               groups.vector = significant.isoform.tissue.selected$groups.vector)
# dev.off()
# 
# absolute.promoter.activity.mean.tissue.cluster2 = as.data.frame(absolute.promoter.activity.mean.tissue[c(major.cluster2, minor.cluster2),]) #
# absolute.promoter.activity.mean.tissue.cluster2$category = c(rep(1,length(major.cluster2)), rep(2,length(minor.cluster2)))
# pheatmap(absolute.promoter.activity.mean.tissue.cluster2[,-ncol(absolute.promoter.activity.mean.tissue.cluster2)], 
#          color = colorRampPalette(c("#00bcd4","white","#FB7AFC"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=10, legend = TRUE,
#          annotation_row=absolute.promoter.activity.mean.tissue.cluster2[, "category", drop=FALSE], 
#          filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.heatmap.down.pdf")))
# 
# 
# ######profileplot and heatmap for different trend isoforms 
# ################################################
# cluster1 = expand.grid(up.cluster, up.cluster.opposite) ##major up regulated
# cluster2 = expand.grid(down.cluster,down.cluster.opposite) ##major down regulated
# build.extract.cluster.name = function(cluster){
#   y = apply(cluster, 1, function(x){getDSPatterns(significant.isoform.tissue.table, x[1], x[2])})
#   return(unlist(y))}
# cluster1.genename = build.extract.cluster.name(cluster1)
# cluster1.genename = cluster1.genename[cluster1.genename %in% significant.isoform.tissue.majoriso.change$L]
# cluster2.genename = build.extract.cluster.name(cluster2)
# cluster2.genename = cluster2.genename[cluster2.genename %in% significant.isoform.tissue.majoriso.change$L]
# 
# ##Major and minor promter name belong to extracted genes
# major.cluster1 = as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,"Major"])
# major.cluster2 = as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,"Major"])
# build.strsplit = function(x){return(unique(unlist(strsplit(x,"_"))))}
# minor.cluster1 = build.strsplit(as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,"Minor"]))
# minor.cluster2 = build.strsplit(as.character(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,"Minor"]))
# write.table(significant.isoform.tissue.table$IsoMajorMinor[cluster1.genename,], file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.mixedup.txt")), sep="\t", col.names = NA, row.names = TRUE, quote = FALSE)
# write.table(significant.isoform.tissue.table$IsoMajorMinor[cluster2.genename,gender,], file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("majoriso.change.mixeddown.txt")), sep="\t", col.names = NA, row.names = TRUE, quote = FALSE)
# 
# ##Plot Major and minor profile in one plot
# load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
# significant.isoform.tissue.selected = significant.isoform.tissue$get2$sig.genes
# significant.isoform.tissue.selected.minor.cluster1 = significant.isoform.tissue.selected$sig.profiles[rownames(significant.isoform.tissue.selected$sig.profiles) %in% minor.cluster1,]
# significant.isoform.tissue.selected.minor.cluster2 = significant.isoform.tissue.selected$sig.profiles[rownames(significant.isoform.tissue.selected$sig.profiles) %in% minor.cluster2,]
# significant.isoform.tissue.selected.major.cluster1 = significant.isoform.tissue.selected$sig.profiles[rownames(significant.isoform.tissue.selected$sig.profiles) %in% major.cluster1,]
# significant.isoform.tissue.selected.major.cluster2 = significant.isoform.tissue.selected$sig.profiles[rownames(significant.isoform.tissue.selected$sig.profiles) %in% major.cluster2,]
# 
# pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.cluster.mixedup.pdf")))
# PlotTwoGroups(data = significant.isoform.tissue.selected.major.cluster1, 
#               data2 = significant.isoform.tissue.selected.minor.cluster1, 
#               x.labels = unique(sample.info$stage),
#               edesign = significant.isoform.tissue.selected$edesign, 
#               groups = significant.isoform.tissue.selected$edesign[, c(3:ncol(significant.isoform.tissue.selected$edesign)),drop=FALSE], 
#               summary.mode = "mean",
#               groups.vector = significant.isoform.tissue.selected$groups.vector)
# dev.off()
# pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.cluster.mixeddown.pdf")))
# PlotTwoGroups(data = significant.isoform.tissue.selected.major.cluster2, 
#               data2 = significant.isoform.tissue.selected.minor.cluster2,
#               x.labels = unique(sample.info$stage),
#               edesign = significant.isoform.tissue.selected$edesign, 
#               groups = significant.isoform.tissue.selected$edesign[, c(3:ncol(significant.isoform.tissue.selected$edesign)),drop=FALSE], 
#               summary.mode = "mean",
#               groups.vector = significant.isoform.tissue.selected$groups.vector)
# dev.off()
# 
# ##Plot Major and minor heatmap
# load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.RData")))
# load(file.path(prewd,"ProcessedData","PromoterObjects","absolute.promoter.activity.mean.RData"))  
# absolute.promoter.activity.mean.tissue = absolute.promoter.activity.mean[,grep(name, colnames(absolute.promoter.activity.mean))]
# absolute.promoter.activity.mean.tissue.cluster1 = as.data.frame(absolute.promoter.activity.mean.tissue[c(major.cluster1, minor.cluster1),]) #
# absolute.promoter.activity.mean.tissue.cluster1$category = c(rep(1,length(major.cluster1)), rep(2,length(minor.cluster1)))
# absolute.promoter.activity.mean.tissue.cluster2 = as.data.frame(absolute.promoter.activity.mean.tissue[c(major.cluster2, minor.cluster2),]) #
# absolute.promoter.activity.mean.tissue.cluster2$category = c(rep(1,length(major.cluster2)), rep(2,length(minor.cluster2)))
# 
# pheatmap(absolute.promoter.activity.mean.tissue.cluster1[,-ncol(absolute.promoter.activity.mean.tissue.cluster1)], 
#          color = colorRampPalette(c("#00bcd4","white","#FB7AFC"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=10, legend = TRUE,
#          annotation_row=absolute.promoter.activity.mean.tissue.cluster1[, "category", drop=FALSE], 
#          filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.heatmap.mixedup.pdf")))
# #annotation_colors = as.list(absolute.promoter.activity.mean.tissue.cluster1[,"color"]),
# pheatmap(absolute.promoter.activity.mean.tissue.cluster2[,-ncol(absolute.promoter.activity.mean.tissue.cluster2)], 
#          color = colorRampPalette(c("#00bcd4","white","#FB7AFC"))(100), border_color = NA, cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE,show_colnames = T,scale = "row", width = 5, height=10, legend = TRUE,
#          annotation_row=absolute.promoter.activity.mean.tissue.cluster2[, "category", drop=FALSE], 
#          filename = file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("majoriso.change.heatmap.mixeddown.pdf")))
# 
# }

# sankey plot
# significant.isoform.tissue.table.extract = as.data.frame(significant.isoform.tissue.table$IsoTable)
# significant.isoform.tissue.table.extract$Cluster.minor = gsub("_[0-9]", "", gsub("0_", "", significant.isoform.tissue.table.extract$Cluster.minor))
# significant.isoform.tissue.table.plot = aggregate(Freq ~ Cluster.Mayor+Cluster.minor,data=significant.isoform.tissue.table.extract,FUN=sum)
# significant.isoform.tissue.table.plot$type = paste("Major",significant.isoform.tissue.table.plot$Cluster.Mayor,"vs_Minor", significant.isoform.tissue.table.plot$Cluster.minor,sep ="_")
# significant.isoform.tissue.table.plot$Freq = as.numeric(significant.isoform.tissue.table.plot$Freq)
# build.up.down.cluster = function(x){
#   if (all(x[1:2] %in% as.character(up.cluster))){y = "Up"} else if (all(x[1:2] %in% as.character(down.cluster))) {y = "Down"} else if (all(x[1:2] %in% c(up.cluster,flat.cluster))) {y = "Up0"} else if (all(x[1:2] %in% c(down.cluster,flat.cluster))) {y = "Down0"} else {y = "Mixed"}
#   return(y)
# }
# significant.isoform.tissue.table.plot$category = factor(apply(significant.isoform.tissue.table.plot,1,build.up.down.cluster), levels = c("Up0", "Up","Down0","Down","Mixed"))
# significant.isoform.tissue.table.plot = significant.isoform.tissue.table.plot[order(significant.isoform.tissue.table.plot$category, significant.isoform.tissue.table.plot$Cluster.Mayor),]


# c(sig.plot[sig.plot$Cluster.Mayor %in% as.character(up.cluster) & sig.plot$Cluster.minor %in% 0, "Freq"],
#   sig.plot[sig.plot$Cluster.Mayor %in% as.character(down.cluster) & sig.plot$Cluster.minor %in% 0, "Freq"],
#   sig.plot[sig.plot$Cluster.Mayor %in% as.character(up.cluster) & sig.plot$Cluster.minor %in% as.character(up.cluster), "Freq"],
#   sig.plot[sig.plot$Cluster.Mayor %in% as.character(down.cluster) & sig.plot$Cluster.minor %in% as.character(down.cluster), "Freq"],
#   sig.plot[sig.plot$Cluster.Mayor %in% as.character(up.cluster) & sig.plot$Cluster.minor %in% as.character(down.cluster), "Freq"],
#   sig.plot[sig.plot$Cluster.Mayor %in% as.character(down.cluster) & sig.plot$Cluster.minor %in% as.character(up.cluster), "Freq"])