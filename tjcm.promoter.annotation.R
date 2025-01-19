##Mouse
prewd="Mammalian2019Promoter/CM2019RNAseqMouse"
# prewd="Mammalian2019Promoter/fastq.snakepipe.file.out"
gtf.file = "/data/repository/organisms/GRCm38_ensembl/Ensembl/release-91/genes.gtf"
species <- 'Mus_musculus'

###First time
# sample.sheet = read.delim(file.path(prewd,"sample.sheet.txt"),header=TRUE, row.names = 1, stringsAsFactors = FALSE)
# sample.sheet$samplename=rownames(sample.sheet)
# if(TRUE){
#  sample.sheet$stage = sapply(strsplit(sample.sheet$samplename,"_"),"[",2)
#  sample.sheet$tissue = sapply(strsplit(sample.sheet$samplename,"_"),"[",1)
#  sample.sheet[sample.sheet$stage=="s0dpb","stage"]="P0"
#  sample.sheet$condition = paste(sample.sheet$tissue, sample.sheet$stage, sep="_")
#  sample.sheet$replicate = unlist(sapply(split(sample.sheet$samplename,sample.sheet$condition),function(x){return(1:length(x))})[unique(sample.sheet$condition)])
#  sample.sheet$newname = paste0(sample.sheet$condition,"_rep", sample.sheet$replicate)
#  sample.sheet$gender = sapply(strsplit(sample.sheet$samplename,"_"),"[",3)
#  sample.sheet$conditiongender = paste(sample.sheet$condition, sample.sheet$gender, sep="_")
#  sample.sheet$replicategender = unlist(sapply(split(sample.sheet$samplename,sample.sheet$conditiongender),function(x){return(1:length(x))})[unique(sample.sheet$conditiongender)])
#  sample.sheet$newnamegender = paste0(sample.sheet$conditiongender,"_rep", sample.sheet$replicategender)
#  rownames(sample.sheet) = sample.sheet$samplename
#  sample.sheet = sample.sheet[order(sample.sheet$tissue, sample.sheet$stage),]
#  sample.sheet=sample.sheet[, -c(1:4)]
#  sample.sheet.male=subset(sample.sheet$conditiongender,grepl("_male",sample.sheet$conditiongender))
#  sample.sheet.female=subset(sample.sheet$conditiongender,grepl("_female",sample.sheet$conditiongender))
#  write.table(sample.sheet, file = file.path(prewd,"sample.sheet.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
#  write.table(sample.sheet.male, file = file.path(prewd,"sample.sheet.male.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
#  write.table(sample.sheet.female, file = file.path(prewd,"sample.sheet.female.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
#  }


sample.info = read.delim(file.path(prewd,"sample.sheet.txt"),header=TRUE, row.names = 1, stringsAsFactors = FALSE)
mc.cores= 20
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

basic.objects=c("basic.objects",ls())

########################################################
####Creating your own promoter annotations using ProActiv
########################################################
if(TRUE){
  dir.create(file.path(prewd,"ProcessedData"))
  dir.create(file.path(prewd,"ProcessedData","AnnotationObjects"))
  
  # make txdb object from gtf file (Gencode v38 used here)
  system(paste("grep -v retained_intron",gtf.file, "|grep -v nonsense_mediated_decay >", c(file.path(prewd,"ProcessedData","AnnotationObjects",'protein_coding.gtf')),sep=" "))
  # db = ensDbFromGtf(gtf = file.path(prewd,"ProcessedData","AnnotationObjects",'protein_coding.gtf'),path=file.path(prewd,"ProcessedData","AnnotationObjects"),organism="Homo_sapiens",version=91,genomeVersion="GRCh38")
  # ensdb = EnsDb(db)
  txdb = makeTxDbFromGFF(file.path(prewd,"ProcessedData","AnnotationObjects",'protein_coding.gtf'))
  saveDb(txdb,file.path(prewd,"ProcessedData","AnnotationObjects",'txdb.sqlite'))
  
  # The species argument to be used for GenomeInfoDb::keepStandardChromosomes
  txdb <- loadDb(file.path(prewd,"ProcessedData","AnnotationObjects",'txdb.sqlite'))
  
  ### Annotation data preparation, delete internalPromoter which overlapped with inner exon inside gene 
  ### Needs to be executed once per annotation. Results can be saved and loaded later for reuse
  promoterAnnotationData <- preparePromoterAnnotation(txdb=txdb, species = species)
  # promoterCoordinates(promoterAnnotationData)[is.na(promoterCoordinates(promoterAnnotationData)$internalPromoter)]$internalPromoter = FALSE
  # promoterCoordinates(promoterAnnotationData) = promoterCoordinates(promoterAnnotationData)[!promoterCoordinates(promoterAnnotationData)$internalPromoter]
  # promoterIdMapping(promoterAnnotationData) = promoterIdMapping(promoterAnnotationData)[promoterIdMapping(promoterAnnotationData)$promoterId %in% promoterCoordinates(promoterAnnotationData)$promoterId,]
  # reducedExonRanges(promoterAnnotationData) = reducedExonRanges(promoterAnnotationData)[reducedExonRanges(promoterAnnotationData)$promoterId %in% promoterCoordinates(promoterAnnotationData)$promoterId]
  # annotatedIntronRanges(promoterAnnotationData) = NULL #annotatedIntronRanges(promoterAnnotationData)[annotatedIntronRanges(promoterAnnotationData)$promoterId %in% promoterCoordinates(promoterAnnotationData)$promoterId]
  
  save(promoterAnnotationData,file=file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
  load(file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
  promoter.anno = promoterAnnotationData@promoterCoordinates
  promoter.anno$name = paste(promoter.anno$promoterId, promoter.anno$geneId, sep="_")
  rtracklayer::export(promoter.anno, con = file.path(prewd,"ProcessedData","AnnotationObjects","promoter.bed"))
  # Retrieve the id mapping between transcripts, TSSs, promoters and genes
  head(promoterIdMapping(promoterAnnotationData))
  
  # Retrieve promoter coordinates
  head(promoterCoordinates(promoterAnnotationData))
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}


################################################
#####Estimate Promoter Activity (STAR alignment) with ProActiv
################################################
if(TRUE){
  method="Proactiv"
  # Preprocessed data is available as part of the package for the human genome (hg19):
  # Available data: proActiv::promoterAnnotationData.gencode.v19
  
  dir.create(file.path(prewd,"ProcessedData","PromoterObjects"))
  dir.create(file.path(prewd,"ProcessedData","PromoterObjects",method))
  
  ### STAR Junction Files Example
  
  # The paths and labels for samples
  starJunctionFiles <- Sys.glob(file.path(prewd,"STAR","*/*.SJ.out.tab"))
  starJunctionFileLabels <- gsub(".SJ.out.tab","",sapply(strsplit(starJunctionFiles,"/"),function(x){return(x[8])}))
  
  
  ## row and column name selected, remove internal promoter and NA promoter
  load(file.path(prewd,"ProcessedData","AnnotationObjects","promoterAnnotationData.RData"))
  promoterCoordinates(promoterAnnotationData)[is.na(promoterCoordinates(promoterAnnotationData)$internalPromoter)]$internalPromoter = TRUE
  promoterCoordinates(promoterAnnotationData) = promoterCoordinates(promoterAnnotationData)[!promoterCoordinates(promoterAnnotationData)$internalPromoter]
  
  # Count the total number of junction reads for each promoter
  promoterCounts.star <- calculatePromoterReadCounts(promoterAnnotationData,
                                                     junctionFilePaths = starJunctionFiles,
                                                     junctionFileLabels =  starJunctionFileLabels,
                                                     junctionType = "star",numberOfCores = mc.cores)
  promoterCounts.star = promoterCounts.star[as.character(promoterCoordinates(promoterAnnotationData)$promoterId),]
  
  ##change counts sample names
  promoterCounts.star = promoterCounts.star[,sample.info$samplename]
  colnames(promoterCounts.star) = sample.info[colnames(promoterCounts.star), "newnamegender"]
  save(promoterCounts.star,file=file.path(prewd,"ProcessedData","PromoterObjects",method,"promoterCounts.star.RData"))
}


########################################################
####DEXseq linux commands
########################################################
##Human
##Exon annotation and count reads with DEXseq and featurecounts
# cd Mammalian2019Promoter/CM2019RNAseqHuman/ProcessedData
# module load HTSeq/0.7.2
# SlurmEasy "python ~/scripts/dexseq_prepare_annotation2.py -f AnnotationObjects/protein_coding.DEXSeq.featurecount.gtf -r yes AnnotationObjects/protein_coding.gtf AnnotationObjects/protein_coding.DEXSeq.gff"
# 
# module load RSeQC/2.6.4
# for i in ../STAR/*bam; do j=`echo $i |sed 's/.bam/.strand.txt/g'`; SlurmEasy "infer_experiment.py -r /data/repository/organisms/GRCh38_ensembl/ensembl/release-91/genes.bed -i $i > $j"; done
# 
# module load subread/1.5.3 && mkdir -p PromoterObjects/Dexseq
# SlurmEasy -t 64 "featureCounts -f -O -s 2 -p -T 64 -F GTF -a AnnotationObjects/protein_coding.DEXSeq.featurecount.gtf -o PromoterObjects/Dexseq/dexseq.counts.txt ../STAR/*bam"

#Mouse
#Exon annotation and count reads with DEXseq and featurecounts
# cd Mammalian2019Promoter/CM2019RNAseqMouse/ProcessedData
# module load HTSeq/0.7.2
# SlurmEasy "python ~/scripts/dexseq_prepare_annotation2.py -f AnnotationObjects/protein_coding.DEXSeq.featurecount.gtf -r yes AnnotationObjects/protein_coding.gtf AnnotationObjects/protein_coding.DEXSeq.gff"
# 
# module load RSeQC/2.6.4
# for i in ../STAR/*bam; do j=`echo $i |sed 's/.bam/.strand.txt/g'`; SlurmEasy "infer_experiment.py -r /data/repository/organisms/GRCm38_ensembl/ensembl/release-91/genes.bed -i $i > $j"; done
# 
# module load subread/1.5.3 && mkdir -p PromoterObjects/Dexseq
# SlurmEasy -t 64 "featureCounts -f -O -s 2 -p -T 64 -F GTF -a AnnotationObjects/protein_coding.DEXSeq.featurecount.gtf -o PromoterObjects/Dexseq/dexseq.counts.txt ../STAR/*bam"
# #SlurmEasy -t 64 "featureCounts -p -B -C -Q 10 --primary -T 8 -s 2 -a ../Annotation/genes.filtered.gtf -o counts.ori.tsv ../STAR/*bam" #only for mouse

##DEX-seq is purely for verification of the method Proativ
########################################################
####Creating your own promoter annotations using DEXseq
########################################################
if(TRUE){
  ##check if the assgined exon name is the same as the original flattened gtf file
  flattenedfile = file.path(prewd, "ProcessedData", "AnnotationObjects", "protein_coding.DEXSeq.featurecount.gtf")
  aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE,header = FALSE)
  colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
  aggregates$strand <- gsub("\\.", "*", aggregates$strand)
  aggregates <- aggregates[which(aggregates$class == "exon"),] # exonic_part
  aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
  aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1",  aggregates$attr)
  # trim the gene_ids to 255 chars in order to match with featurecounts
  longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
  warning(paste0(longIDs, " aggregate geneIDs were found truncated in featureCounts output"), call. = FALSE)
  aggregates$gene_id <- substr(aggregates$gene_id,1,255)
  
  exonids <- gsub(".*exon_number\\s(\\S+).*", "\\1",aggregates$attr) # exonic_part_number
  exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start,end = aggregates$end), strand = aggregates$strand)
  names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")
  transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
  transcripts <- strsplit(transcripts, "\\+")
  names(transcripts) <- names(exoninfo)
  # if (!all(rownames(dcounts) %in% names(exoninfo))) {
  #   stop("Count files do not correspond to the flattened annotation file")
  # }
  # matching <- match(rownames(dcounts), names(exoninfo))
  # stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
  # stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
  rtracklayer::export(exoninfo, con = file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exon.bed"))
  save(exoninfo, file = file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.RData"))
  save(transcripts, file = file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.transcripts.RData"))
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

#########################################################
###DEXSeq Exon annotation correpond promoter annotation
#########################################################
if(TRUE){
  load(file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.RData"))
  load(file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
  promoter.annotation = promoterAnnotationData@promoterCoordinates
  seqlevelsStyle(promoter.annotation) = "Ensembl"
  
  dexseq.exoninfo.promoter.overlap = findOverlaps(promoter.annotation, exoninfo)
  dexseq.exoninfo.promoter.correspond = as.data.frame(mcols(promoter.annotation[queryHits(dexseq.exoninfo.promoter.overlap)])[,c("promoterId","geneId","internalPromoter")], stringsAsFactors = FALSE)
  dexseq.exoninfo.promoter.correspond$dexseq = names(exoninfo[subjectHits(dexseq.exoninfo.promoter.overlap)])
  save(dexseq.exoninfo.promoter.correspond, file = file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.promoter.correspond.RData"))
  write.table(dexseq.exoninfo.promoter.correspond, file = file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.promoter.correspond.txt"), row.names = FALSE, quote = FALSE, sep="\t")
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

################################################
#####Estimate Promoter Activity (STAR alignment) with DEXseq
################################################
if(TRUE){
  method="Dexseq"
  dir.create(file.path(prewd,"ProcessedData","PromoterObjects"))
  dir.create(file.path(prewd,"ProcessedData","PromoterObjects",method))
  
  ##dexseq counts matrix from feature-count result
  counts.table = read.delim(file = file.path(prewd, "ProcessedData", "PromoterObjects", method, "dexseq.counts.txt"),comment.char = "#")
  colnames(counts.table) = gsub(".bam", "", gsub("[.]*STAR[.]", "", colnames(counts.table)))
  
  ##change rownames
  dcounts = counts.table
  id <- as.character(dcounts[,1])
  n <- id
  split(n,id) <- lapply(split(n ,id), seq_along )
  rownames(dcounts) <- sprintf("%s%s%03.f",id,":E",as.numeric(n))
  dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginnning of gene name 
  dcounts <- dcounts[,7:ncol(dcounts)]
  
  ##change colnames
  dcounts = dcounts[,sample.info$samplename]
  colnames(dcounts) = sample.info[colnames(dcounts), "newnamegender"]

  save(dcounts, file= file.path(prewd, "ProcessedData", "PromoterObjects", method, "dexseq.counts.dcounts.RData"))
  
  ##extract only promoters
  load(file.path(prewd,"ProcessedData","PromoterObjects", "Proactiv",'promoterCounts.star.RData'))
  load(file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.promoter.correspond.RData"))
  dcounts = dcounts[,colnames(promoterCounts.star)]
  dcounts.promoter = merge(dexseq.exoninfo.promoter.correspond[,c("promoterId","dexseq")], dcounts, by.x ="dexseq", by.y = 0)
  rownames(dcounts.promoter) = dcounts.promoter$promoterId
  dcounts.promoter = dcounts.promoter[,-1:-2]
  
  ## remove internal promoter and NA promoter
  load(file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
  promoterCoordinates(promoterAnnotationData)[is.na(promoterCoordinates(promoterAnnotationData)$internalPromoter)]$internalPromoter = TRUE
  promoterCoordinates(promoterAnnotationData) = promoterCoordinates(promoterAnnotationData)[!promoterCoordinates(promoterAnnotationData)$internalPromoter]
  promoterCounts.star = dcounts.promoter[promoterCoordinates(promoterAnnotationData)$promoterId,]
  save(promoterCounts.star,file=file.path(prewd,"ProcessedData","PromoterObjects", method,'promoterCounts.star.RData'))
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

#########################################################
###DEXSeq Exon counts versus promoter activity
#########################################################
if(TRUE){
  # load(file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.promoter.correspond.RData"))
  proactiv = get(load(file.path(prewd, "ProcessedData", "PromoterObjects", "Proactiv", "promoterCounts.star.RData")))
  dexseq = get(load(file.path(prewd,"ProcessedData","PromoterObjects", "Dexseq","promoterCounts.star.RData")))
  
  tmp1 = proactiv[,1,drop=FALSE]
  tmp2 = dexseq[,1,drop=FALSE]
  tmp1 = tmp1[tmp1[,1]>=10,,drop=FALSE]
  tmp2 = tmp2[rownames(tmp1),,drop=FALSE]
  plot(log2(tmp1[,1]+1),log2(tmp2[,1]+1))
  cor(log2(tmp1[,1]+1),log2(tmp2[,1]+1))
  

  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}



################################################
###counts normalization size factor
################################################
if(TRUE){
  dir.create(file.path(prewd,"bamCoverage"))
  ## load dexseq counts of promoters
  dcounts.promoter = get(load(file.path(prewd,"ProcessedData","PromoterObjects","Dexseq",'promoterCounts.star.RData')))  #dcounts[dexseq.exoninfo.promoter.correspond[dexseq.exoninfo.promoter.correspond$promoterId %in% rownames(promoterCounts.star),"dexseq"],]
  
  ## load junction counts of promoters
  load(file.path(prewd,"ProcessedData","PromoterObjects","Proactiv","promoterCounts.star.RData"))
  
  ## load DEXSeq counts of exons
  load(file.path(prewd, "ProcessedData", "AnnotationObjects", "dexseq.exoninfo.promoter.correspond.RData"))
  load(file.path(prewd, "ProcessedData", "PromoterObjects","Dexseq", "dexseq.counts.dcounts.RData"))
  dcounts = dcounts[,colnames(promoterCounts.star)]
  
  ## load feature counts of genes
  feature.counts = read.delim(file.path(prewd,"featureCounts","counts.tsv"),header=TRUE, row.names = 1, stringsAsFactors = FALSE)
  colnames(feature.counts) = sample.info[colnames(feature.counts), "newnamegender"]
  feature.counts = feature.counts[,colnames(promoterCounts.star)]
  
  ##DEseq2
  build.deseq.sf = function(counts){
    colData <- data.frame(sampleLabels = colnames(counts))
    rownames(colData) <- colnames(counts)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[which(!is.na(counts[, 1])),], colData = colData, design = ~1)
    dds <- DESeq2::estimateSizeFactors(dds)
    sf.deseq2 = 1/sizeFactors(dds)
    return(sf.deseq2)
  }
  junction.promoter.deseq2.sf = as.data.frame(build.deseq.sf(promoterCounts.star))
  dexseq.promoter.deseq2.sf = as.data.frame(build.deseq.sf(dcounts.promoter))
  dexseq.exons.deseq2.sf = as.data.frame(build.deseq.sf(dcounts))
  feature.counts.deseq2.sf = as.data.frame(build.deseq.sf(feature.counts))
  cor(dexseq.exons.deseq2.sf, junction.promoter.deseq2.sf)
  cor(dexseq.exons.deseq2.sf, dexseq.promoter.deseq2.sf)
  cor(dexseq.exons.deseq2.sf, feature.counts.deseq2.sf)
  # plot(sf.deseq2, sf.tmm, xlim=c(0,2), ylim=c(0,2))
  
  ##save DESeq2 all exons and all genes normalization
  rownames(sample.info) = sample.info$newnamegender
  dexseq.promoter.deseq2.sf$oriname = sample.info[rownames(dexseq.promoter.deseq2.sf), "samplename"]
  dexseq.exons.deseq2.sf$oriname = sample.info[rownames(dexseq.exons.deseq2.sf), "samplename"]
  feature.counts.deseq2.sf$oriname = sample.info[rownames(feature.counts.deseq2.sf), "samplename"]
  
  # write.table(dexseq.promoter.deseq2.sf,file=file.path(prewd,"bamCoverage/deseq2.dexseq.promoter.txt"), sep =":", quote = FALSE, row.names = TRUE, col.names = NA)
  # write.table(dexseq.exons.deseq2.sf,file=file.path(prewd,"bamCoverage/deseq2.dexseq.counts.txt"), sep =":", quote = FALSE, row.names = TRUE, col.names = NA)
  # write.table(feature.counts.deseq2.sf,file=file.path(prewd,"bamCoverage/deseq2.feature.counts.txt"), sep =":", quote = FALSE, row.names = TRUE, col.names = NA)
  
  ##edgeR
  build.edger.sf = function(counts){
    cds = DGEList(counts=as.matrix(counts[which(!is.na(counts[, 1])),]),group=factor(colnames(counts))) #1:ncol(counts)
    cds = calcNormFactors(cds, method=c("TMM"))
    sf.tmm = cds$samples$norm.factors
    names(sf.tmm) = rownames(cds$samples)
    return(1/sf.tmm)
  }
  junction.promoter.edger.sf = as.data.frame(build.edger.sf(promoterCounts.star))
  dexseq.promoter.edger.sf = as.data.frame(build.edger.sf(dcounts.promoter))
  dexseq.exons.edger.sf = as.data.frame(build.edger.sf(dcounts))
  feature.counts.edger.sf = as.data.frame(build.edger.sf(feature.counts))
  cor(dexseq.exons.edger.sf, junction.promoter.edger.sf)
  cor(dexseq.exons.edger.sf, dexseq.promoter.edger.sf)
  cor(dexseq.exons.edger.sf, feature.counts.edger.sf)
  
  ##save edgeR all exons and all genes normalization
  rownames(sample.info) = sample.info$newnamegender
  dexseq.promoter.edger.sf$oriname = sample.info[rownames(dexseq.promoter.edger.sf), "samplename"]
  dexseq.exons.edger.sf$oriname = sample.info[rownames(dexseq.exons.edger.sf), "samplename"]
  feature.counts.edger.sf$oriname = sample.info[rownames(feature.counts.edger.sf), "samplename"]
  
  # write.table(dexseq.promoter.edger.sf,file=file.path(prewd,"bamCoverage/edger.dexseq.promoter.txt"), sep =":", quote = FALSE, row.names = TRUE, col.names = NA)
  # write.table(dexseq.exons.edger.sf,file=file.path(prewd,"bamCoverage/edger.dexseq.counts.txt"), sep =":", quote = FALSE, row.names = TRUE, col.names = NA)
  # write.table(feature.counts.edger.sf,file=file.path(prewd,"bamCoverage/edger.feature.counts.txt"), sep =":", quote = FALSE, row.names = TRUE, col.names = NA)
  
  # plot(sf.deseq2, sf.tmm, xlim=c(0,2), ylim=c(0,2))
  
  ##DEXSeq
  # sample.sheet = data.frame(row.names = colnames(dcounts), condition = gsub("_.*", "", colnames(dcounts)), stringsAsFactors = FALSE)
  # dxd = DEXSeqDataSet(countData = dcounts, sampleData = sample.sheet, design = ~sample + exon, featureID = sapply(strsplit(rownames(dcounts), ":"), "[[", 2), groupID = sapply(strsplit(rownames(dcounts), ":"), "[[", 1) ) #,
  # #, featureRanges = exoninfo[matching], transcripts = transcripts[matching]
  # # dxd = dxd[rownames(dxd) %in% dexseq.exoninfo.promoter.correspond$dexseq,]  #take only promoters
  # dxd = dxd[!apply(featureCounts(dxd), 1, function(x){return(all(x==0))}),]  #remove exons with all 0s
  # dxd = estimateSizeFactors( dxd )
  # dexseq.exons.dexseq.sf = 1/sizeFactors(dxd)[1:(length(sizeFactors(dxd))/2)] 
  ##Turns out it is exactly the same as dexseq.exons.dexseq.sf
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}




