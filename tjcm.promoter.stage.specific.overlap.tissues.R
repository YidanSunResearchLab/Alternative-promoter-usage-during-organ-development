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
#library(plotly)
#library(qdapTools)
#library(UpSetR)

method="Proactiv"
# method="Dexseq"
norm.method = "edger"
# norm.method = "deseq2"
gender = "Male"
# gender = "Female"
# gender = "ale"
category.order = c("Down_Up","Up_Down","Flat_Down","Down_Flat","Down_Down","Flat_Up","Up_Flat","Up_Up")

source("~/scripts/tjcm.promoter.stage.specific.function.R")

if(gender == "Male"){
  tissue.order = c("Brain", "Cerebellum", "Heart", "Kidney" , "Liver", "Testis","Ovary")
} else if (gender == "Female") {
  tissue.order = c("Ovary") #"Brain", "Cerebellum", "Heart", "Kidney" , "Liver", 
} else {
  tissue.order = c("Brain", "Cerebellum", "Heart", "Kidney" , "Liver", "Ovary", "Testis")
}
mc.cores= 20
flat.cluster = c(0)
up.cluster = c(2) #2,3,6
down.cluster = c(1) #1,4,5
theme_classic.add = theme_classic()+ theme(legend.position = "none",strip.background = element_blank(),strip.text = element_blank(), plot.title = element_text(hjust = 0.5, family = "Helvetica"),axis.text = element_text(colour = "black", family = "Helvetica"), axis.title=element_text(colour = "black", family = "Helvetica"), text = element_text(colour = "black", family = "Helvetica"))
theme_bw.add = theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, family = "Helvetica"),axis.text = element_text(colour = "black", family = "Helvetica"), axis.title=element_text(colour = "black", family = "Helvetica"), text = element_text(colour = "black", family = "Helvetica"))

basic.objects=c("basic.objects",ls(),"name")

#############################################################
##This script is mainly for compare different tissues
##Check the overlap between different tissues
##Part 1: promoter number statistics in Up/Down cluster in each tissue
##Part 2: up cluster overlap between all tissues, so for down and flat cluster
##Part 3: table categories of overlapped genes and tissue specific genes
##Part 3: GO annotation of overlapped genes and tissue specific genes
############Three important table for all tissues overlap #####################
##cluster category (up,down,flat) info: all.tissues.promoter.cluster.category.list
##table category info: all.tissues.promoter.table.category.list ("Down_Down","Down_Flat","Flat_Down","Down_Up","Up_Down","Flat_Up","Up_Flat","Up_Up")
##overlap between tissues: all.tissues.promoter.up/down/flat.cluster.overlap
############Two important table for each tissueafterwards#####################
##in each tissue, find DE promoters common with others: all.tissues.promoter.table.category.common (at least two tissues have this)
##in each tissue, find DE promoters specific to this tissue: all.tissues.promoter.table.category.specific

############################################################
##Prepare cluster info and table category info from all tissues 
############################################################
if (TRUE){
    
  dir.create(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration"))
  dir.create(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender))
  dir.create(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration"))
  dir.create(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender))
  
  all.tissues.promoter.table.category.list = list()
  all.tissues.promoter.cluster.category.list = list()
  
  for (name in rev(tissue.order)){
    if(name=="Ovary"){gender="Female"} else {gender="Male"}
    message(name)
    # name="Testis"
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster0.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.category.RData")))
    
    ##The cluster name is flipped for Kidney, usually 2 is up. However 1 is up for Kidney
    if(name=="Kidney"){ 
      #flipping_map <- c(
      #  "Up_Up" = "Down_Down",
      #  "Down_Down" = "Up_Up",
      #  "Flat_Up" = "Flat_Down",
      #  "Flat_Down" = "Flat_Up",
      #  "Up_Flat" = "Down_Flat",
      #  "Down_Flat" = "Up_Flat",
      #  "Up_Down" = "Down_Up",
      #  "Down_Up" = "Up_Down"
      #)
      
      flipped_number = c("1"=2, "2"=1)
      
      # Apply the flipping map
      #significant.isoform.tissue.table.category$category = flipping_map[significant.isoform.tissue.table.category$category]
      significant.isoform.tissue.see.cluster$cluster= as.numeric(flipped_number[significant.isoform.tissue.see.cluster$cluster])
    }
    
    ##merge cluster info and table category info from different tissues into one list
    all.tissues.promoter.table.category.list[[name]] = significant.isoform.tissue.table.category
    all.tissues.promoter.cluster.category.list[[name]] = rbind(significant.isoform.tissue.see.cluster, significant.isoform.tissue.see.cluster0)
    print(table( all.tissues.promoter.table.category.list[[name]]$category))
    print(table( all.tissues.promoter.cluster.category.list[[name]]$cluster))
    
  }
  
  
  save(all.tissues.promoter.table.category.list, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.table.category.list.RData"))
  save(all.tissues.promoter.cluster.category.list, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.cluster.category.list.RData"))
  
  
  ##Save all.tissues.promoter.table.category.list and all.tissues.promoter.cluster.category.list which stores all the significant DDPs in each tissue
  save(all.tissues.promoter.table.category.list, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration","all.tissues.promoter.table.category.list.RData"))
  save(all.tissues.promoter.cluster.category.list, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration","all.tissues.promoter.cluster.category.list.RData"))
  
  all.tissues.promoter.cluster.category.list.df <- do.call(rbind, all.tissues.promoter.cluster.category.list)
  all.tissues.promoter.cluster.category.list.df$tissue <- sub("^([A-Za-z]+)\\..*", "\\1", rownames(all.tissues.promoter.cluster.category.list.df))
  all.tissues.promoter.cluster.category.list.df$GeneID <- sub(".*[.]", "", rownames(all.tissues.promoter.cluster.category.list.df))
  write.table(all.tissues.promoter.cluster.category.list.df, file = file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",paste0("all.tissues.ddp.txt")), quote = FALSE, col.names = NA, sep="\t")
  
  all.tissues.promoter.table.category.list.df <- do.call(rbind, all.tissues.promoter.table.category.list)
  all.tissues.promoter.table.category.list.df$tissue <- sub("^([A-Za-z]+)\\..*", "\\1", rownames(all.tissues.promoter.table.category.list.df))
  all.tissues.promoter.table.category.list.df$GeneID <- sub(".*[.]", "", rownames(all.tissues.promoter.table.category.list.df))
  write.table(all.tissues.promoter.table.category.list.df, file = file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",paste0("all.tissues.ddp.category.txt")), quote = FALSE, col.names = NA, sep="\t")
  #Supplementary table S2 in the paper
  ## cd Mammalian2019Promoter/CM2019RNAseqMouse/Output/StageSpecific/Proactiv/edger/Integration/Male/ && cat *ddp.txt > TableS1.txt
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
  
}

##protein isoform analysis
if(TRUE){
  gtf.file = import("/data/repository/organisms/GRCm38_ensembl/Ensembl/release-91/genes.gtf")
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityTable.RData"))
  load(file.path(prewd,"ProcessedData","AnnotationObjects",'promoterAnnotationData.RData'))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration","Male","all.tissues.promoter.table.category.list.RData"))
  
  
  transcript.biotype = unique(as.data.frame(mcols(gtf.file)[,c("transcript_id","transcript_biotype" )]))
  promoter.transcript=promoterIdMapping(promoterAnnotationData)
  promoter.transcript$promoterId = paste0("prmtr.",promoter.transcript$promoterId)
  promoter.transcript.merged <- merge( promoter.transcript,  transcript.biotype,  by.x = "transcriptName", by.y = "transcript_id", all.y = FALSE, all.x = TRUE )
 
  for (name in names(all.tissues.promoter.table.category.list)){
    all.tissues.promoter.table.category.list.selected=all.tissues.promoter.table.category.list[[name]]
    print(name)
    # Load necessary libraries
    library(dplyr)
    
    # Iterate over each category in the first dataframe
    categories <- unique(all.tissues.promoter.table.category.list.selected$category)
    
    # Initialize an empty list to store results for plotting
    biotype_counts <- list()
    
    for (cat in categories) {
      # Filter rows for the current category
      subset_data <- all.tissues.promoter.table.category.list.selected %>% filter(category == cat)
      
      # Extract Minor column, split by "_", and unlist into individual promoters
      promoters <- unlist(strsplit(subset_data$Minor, "_"))
      
      # Match the promoters to promoterId in the second dataframe
      matching_rows <- promoter.transcript.merged %>%  filter(promoterId %in% promoters)
      
      # Count occurrences of transcript_biotype
      biotype_count <- as.data.frame(table(matching_rows$transcript_biotype))
      biotype_count$category <- ifelse(
        biotype_count$Var1 == "protein_coding",
        "protein_coding",
        "non_coding"
      )
      
      # Sum the counts for each category
      biotype_count <- biotype_count %>%
        group_by(category) %>%
        summarise(total_count = sum(Freq), .groups = "drop")
      
      
      # Append to results list
      biotype_counts[[cat]] <- as.data.frame(biotype_count)
    }
    
    # Combine all counts into a single data frame
    biotype_counts_df <- do.call(rbind, lapply(names(biotype_counts), function(name) transform(biotype_counts[[name]], category_name = name)))
    biotype_counts_df$category_name = factor(biotype_counts_df$category_name, levels=c("Up_Up", "Up_Flat", "Flat_Up", "Down_Down", "Down_Flat","Flat_Down", "Up_Down", "Down_Up"))
    
    # Plot the results using ggplot2
    p = ggplot(biotype_counts_df, aes(x = category_name, y = total_count, fill = category)) +
      geom_bar(stat = "identity", position = "fill", alpha=0.7,colour="white" ) +
      labs(title = name,
           x = "",
           y = "Proportion",
           fill = "Category") +
      scale_fill_manual(values = c("non_coding" = "darkgrey", "protein_coding" = "#FF0075")) +
      theme_classic() +
      theme(
        panel.grid = element_blank(),              # Remove grids
        axis.text = element_text(color = "black"), # Make axis text black
        axis.title = element_text(color = "black"),# Make axis titles black
        plot.title = element_text(color = "black", hjust = 0.5), # Center and make title black
        legend.title = element_text(color = "black"), # Make legend title black
        legend.text = element_text(color = "black")  # Make legend text black
      )
    ggsave(plot=p,file=file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.promoter.proteincoding.",name,".pdf")), width=8, height=2)
    
  }
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

############################################################################
######promoter number statistics in Up/Down cluster in each tissue
############################################################################
if(TRUE){
  # name = "Brain" "Testis"
  
  ##load the cluster.category and table.category info
  #load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration","Female","all.tissues.promoter.cluster.category.list.RData"))
  #all.tissues.promoter.cluster.category.list.female = all.tissues.promoter.cluster.category.list
  #load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration","Female","all.tissues.promoter.table.category.list.RData"))
  #all.tissues.promoter.table.category.female = all.tissues.promoter.table.category.list
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration","Male","all.tissues.promoter.cluster.category.list.RData"))
  #all.tissues.promoter.cluster.category.list[["Ovary"]] = all.tissues.promoter.cluster.category.list.female[["Ovary"]]
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration","Male","all.tissues.promoter.table.category.list.RData"))
  #all.tissues.promoter.table.category.list[["Ovary"]] = all.tissues.promoter.table.category.female[["Ovary"]]
  tissue.number.order = c( "Heart","Ovary","Kidney","Cerebellum","Liver","Brain", "Testis","Ovary")
  
  ##number statistics for up down cluster in all tissues
  #up.all = c(as.character(tmp[tmp$category %in% c("Up_Up","Up_Down", "Up_Flat"),"Major"]),unlist(apply(tmp[tmp$category %in% c("Down_Up","Flat_Up","Up_Up"),"Minor",drop=F],1,function(x){strsplit(x,"_")})) ) ##use another method calculating how many down and up promoters
  #sdown.all = c(as.character(tmp[tmp$category %in% c("Down_Down","Down_Up", "Down_Flat"),"Major"]),unlist(apply(tmp[tmp$category %in% c("Up_Down","Flat_Down","Down_Down"),"Minor",drop=F],1,function(x){strsplit(x,"_")})) )  ##use another method calculating how many down and up promoters
  for (name in names(all.tissues.promoter.table.category.list)){
    tmp=all.tissues.promoter.table.category.list[[name]]
    print(name)
    up.all = c(major_up_up=length(tmp[tmp$category=="Up_Up","Major"]), major_up_flat=length(tmp[tmp$category=="Up_Flat","Major"]), major_up_down=length(tmp[tmp$category=="Up_Down","Major"]),
    minor_up_up=length(unlist(apply(tmp[tmp$category=="Up_Up","Minor",drop=F],1,function(x){strsplit(x,"_")}))), minor_flat_up=length(unlist(apply(tmp[tmp$category=="Flat_Up","Minor",drop=F],1,function(x){strsplit(x,"_")}))), minor_down_up=length(unlist(apply(tmp[tmp$category=="Down_Up","Minor",drop=F],1,function(x){strsplit(x,"_")})))
    )
    down.all = c(major_down_down=length(tmp[tmp$category=="Down_Down","Major"]), major_down_flat=length(tmp[tmp$category=="Down_Flat","Major"]), major_down_up=length(tmp[tmp$category=="Down_Up","Major"]),
               minor_down_down=length(unlist(apply(tmp[tmp$category=="Down_Down","Minor",drop=F],1,function(x){strsplit(x,"_")}))), minor_flat_down=length(unlist(apply(tmp[tmp$category=="Flat_Down","Minor",drop=F],1,function(x){strsplit(x,"_")}))), minor_up_down=length(unlist(apply(tmp[tmp$category=="Up_Down","Minor",drop=F],1,function(x){strsplit(x,"_")})))
    )
    print(up.all)
    print(down.all)
  }
  all.tissue.cluster.category = sapply(all.tissues.promoter.cluster.category.list, function(x){return(table(x$cluster))})
  cluster.name = setNames(c("Flat","Up","Down"),c(flat.cluster, up.cluster, down.cluster))
  rownames(all.tissue.cluster.category) = cluster.name[rownames(all.tissue.cluster.category)]
  all.tissue.cluster.category = all.tissue.cluster.category[-1,]
  all.tissue.cluster.category.plot = melt(all.tissue.cluster.category)
  all.tissue.cluster.category.plot$Var1 = factor(all.tissue.cluster.category.plot$Var1, levels=c("Up","Down"))
  all.tissue.cluster.category.plot$Var2 = factor(as.character(all.tissue.cluster.category.plot$Var2), levels=(levels(all.tissue.cluster.category.plot$Var2)))
  colnames(all.tissue.cluster.category.plot) = c("Type","tissue","value")
  
  p1=ggplot(all.tissue.cluster.category.plot,aes(fill=Type, x=tissue,y=value)) +
    geom_bar(stat="identity", alpha=0.7,width=0.7) + #,position="fill"
    geom_text(aes(label=value),position = position_stack(reverse = FALSE, vjust=0.5), size=3.5) +
    scale_fill_manual(values=c("#f76a8c", "#77d8d8")) + #guides(fill=FALSE) + #FDBCFD
    xlab("") +
    ylab("Promoter Number") +
    coord_flip() +
    theme_classic.add + theme(legend.position = c(0.9, 0.8))
  print(p1)
  ggsave(p1,file=file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.cluster.category.pdf")), width=4, height=4)
  
  ##number statistics for up down cluster  in all tissues
  table.category.statistics = sapply(all.tissues.promoter.table.category.list, function(x){return(table(x$category))})
  table.category.statistics = cbind(table.category.statistics, others=rep(0, nrow(table.category.statistics)))
  table.category.statistics.plot = melt(table.category.statistics)
  table.category.statistics.plot$Var1 = factor(table.category.statistics.plot$Var1, levels = rev(category.order))
  table.category.statistics.plot$Var2 = factor(table.category.statistics.plot$Var2, levels = rev(levels(table.category.statistics.plot$Var2)))#c("others", (colnames(table.category.statistics)[order(colSums(table.category.statistics[c(5,7:8),])/colSums(table.category.statistics))])[-ncol(table.category.statistics)]) )
  colnames(table.category.statistics.plot) = c("Type","tissue","value")
  
  p2=ggplot(table.category.statistics.plot,aes(fill=Type, x=tissue,y=value)) +
    geom_bar(stat="identity",position="fill", alpha=0.8,width=0.7) + #,position="fill"
    geom_text(aes(label=value),position = position_fill(reverse = FALSE, vjust=0.5), size=2) +
    scale_fill_manual(values=c( "#f32053", "#f76a8c", "#fcc7d4", rev(c("#ffce8f","#ff9206", "#b4e9e9","#77d8d8","#3ac7c7" ))) ) + #guides(fill=FALSE) +"#deff8b","#8cba51","#216353" 
    scale_y_continuous(labels = scales::percent) +
    xlab("") +
    ylab("") +
    coord_polar(theta="y") +
    theme_classic.add + theme(legend.position = "right", axis.text.x=element_blank())
  print(p2)
  ggsave(p2,file=file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.table.category.pdf")), width=5, height=5)
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

##############################################################
######up cluster overlap between all tissues, also for down and flat cluster (UpSetR plot)
##############################################################
if(TRUE){
  # name = "Brain" "Testis"
  
  ##load the cluster.category and table.category info
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.cluster.category.list.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.table.category.list.RData"))
  
  ##selected cluster
  build.cluster.plot = function(selected.cluster=up.cluster){
    cluster.selected.list = lapply(all.tissues.promoter.cluster.category.list, function(x){return(rownames(x[x$cluster==selected.cluster,,drop=FALSE]))} )
    cluster.selected.list.binary = t(mtabulate(cluster.selected.list)!=0)
    cluster.selected.list.plot = as.data.frame(apply(cluster.selected.list.binary, 2, as.numeric))
    rownames(cluster.selected.list.plot) = rownames(cluster.selected.list.binary)
    return(cluster.selected.list.plot)
  }   
  
  ##up cluster
  cluster.name= "up"
  all.tissues.promoter.up.cluster.overlap = build.cluster.plot(selected.cluster=up.cluster)
  p1 = upset(all.tissues.promoter.up.cluster.overlap, sets.bar.color = c("#fa7d09","#0779e4","#56B4E9","#fc5185","#ffde7d","#1fab89"),order.by = "freq", 
             main.bar.color = c("#fa7d09","#ffde7d","#0779e4","#1fab89","#fc5185","#56B4E9","black", rep("darkgrey",56)),
             nsets=ncol(all.tissues.promoter.up.cluster.overlap) ,keep.order = T, empty.intersections = "on")
  pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.overlap.", cluster.name,".pdf")), width=3, height=3)
  print(p1)
  dev.off()
  save(all.tissues.promoter.up.cluster.overlap, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.up.cluster.overlap.RData"))
  
  ##down cluster
  cluster.name= "down"
  all.tissues.promoter.down.cluster.overlap = build.cluster.plot(selected.cluster=down.cluster)
  p2 = upset(all.tissues.promoter.down.cluster.overlap, sets.bar.color = c("#fa7d09","#1fab89","#0779e4","#ffde7d","#56B4E9","#fc5185"),order.by = "freq", 
             main.bar.color = c("#fa7d09","#1fab89","#ffde7d","#0779e4","darkgrey", "#56B4E9","#fc5185","black", rep("darkgrey",55)),
             nsets=ncol(all.tissues.promoter.down.cluster.overlap) ,keep.order = T, empty.intersections = "on")
  pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.overlap.", cluster.name,".pdf")), width=3, height=3)
  print(p2)
  dev.off()
  save(all.tissues.promoter.down.cluster.overlap, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.down.cluster.overlap.RData"))
  
  ##flat cluster
  cluster.name= "flat"
  all.tissues.promoter.flat.cluster.overlap = build.cluster.plot(selected.cluster=flat.cluster)
  p3 = upset(all.tissues.promoter.flat.cluster.overlap, sets.bar.color = c("#fa7d09","#0779e4","#56B4E9","#fc5185","#ffde7d","#1fab89"),order.by = "freq", 
             main.bar.color = c("#fa7d09","#ffde7d","#0779e4","#1fab89","#fc5185","#56B4E9","black", rep("darkgrey",56)),
             nsets=ncol(all.tissues.promoter.flat.cluster.overlap) ,keep.order = T, empty.intersections = "on")
  pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.overlap.", cluster.name,".pdf")), width=3, height=3)
  print(p3)
  dev.off()
  save(all.tissues.promoter.flat.cluster.overlap, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.flat.cluster.overlap.RData"))
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

################################################
##table categories of overlapped genes and tissue specific genes
##Brain and Cerellum has special treatment
################################################
for (name in tissue.order){
  # name = "Testis" #"Brain" 
  if(name=="Ovary"){gender="Female"} else {gender="Male"}
  
  ##load the table.category and up/down/flat.cluster.overlap info
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.cluster.category.list.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.table.category.list.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.up.cluster.overlap.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.down.cluster.overlap.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.flat.cluster.overlap.RData"))
  
  all.tissues.promoter.cluster.category.selected = all.tissues.promoter.cluster.category.list[[name]]
  all.tissues.promoter.table.category.selected = all.tissues.promoter.table.category.list[[name]]
  table.category.statistics = data.frame(table(all.tissues.promoter.table.category.selected$category))
  
  ##common with other tissues
  all.tissues.promoter.up.cluster.common = all.tissues.promoter.up.cluster.overlap[all.tissues.promoter.up.cluster.overlap[,name]==1 & rowSums(all.tissues.promoter.up.cluster.overlap) > 1, ]
  if(name == "Brain"){all.tissues.promoter.up.cluster.2common = all.tissues.promoter.up.cluster.common[all.tissues.promoter.up.cluster.common$Cerebellum==1,]; all.tissues.promoter.up.cluster.common = all.tissues.promoter.up.cluster.common[all.tissues.promoter.up.cluster.common$Cerebellum!=1,]}
  if(name == "Cerebellum"){all.tissues.promoter.up.cluster.2common = all.tissues.promoter.up.cluster.common[all.tissues.promoter.up.cluster.common$Brain==1,]; all.tissues.promoter.up.cluster.common = all.tissues.promoter.up.cluster.common[all.tissues.promoter.up.cluster.common$Brain!=1,]}
  all.tissues.promoter.down.cluster.common = all.tissues.promoter.down.cluster.overlap[all.tissues.promoter.down.cluster.overlap[,name]==1 & rowSums(all.tissues.promoter.down.cluster.overlap) > 1, ]
  if(name == "Brain"){all.tissues.promoter.down.cluster.2common = all.tissues.promoter.down.cluster.common[all.tissues.promoter.down.cluster.common$Cerebellum==1,]; all.tissues.promoter.down.cluster.common = all.tissues.promoter.down.cluster.common[all.tissues.promoter.down.cluster.common$Cerebellum!=1,]}
  if(name == "Cerebellum"){all.tissues.promoter.down.cluster.2common = all.tissues.promoter.down.cluster.common[all.tissues.promoter.down.cluster.common$Brain==1,]; all.tissues.promoter.down.cluster.common = all.tissues.promoter.down.cluster.common[all.tissues.promoter.down.cluster.common$Brain!=1,]}
  
  all.tissues.promoter.cluster.common = c(rownames(all.tissues.promoter.up.cluster.common), rownames(all.tissues.promoter.down.cluster.common))
  all.tissues.promoter.cluster.common.geneid = unique(all.tissues.promoter.cluster.category.selected[rownames(all.tissues.promoter.cluster.category.selected) %in% all.tissues.promoter.cluster.common, "geneid"])
  all.tissues.promoter.table.category.common = all.tissues.promoter.table.category.selected[rownames(all.tissues.promoter.table.category.selected) %in% all.tissues.promoter.cluster.common.geneid, ]
  all.tissues.promoter.table.category.common.plot = data.frame(table(all.tissues.promoter.table.category.common$category))
  save(all.tissues.promoter.table.category.common, file=file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.common.RData")))
  
  all.tissues.promoter.table.category.common.plot = merge(data.frame(table(all.tissues.promoter.table.category.common$category)),table.category.statistics, by.x="Var1", by.y="Var1")
  all.tissues.promoter.table.category.common.plot$rawvalue = all.tissues.promoter.table.category.common.plot$Freq.x
  all.tissues.promoter.table.category.common.plot$percentage = round(all.tissues.promoter.table.category.common.plot$Freq.x/all.tissues.promoter.table.category.common.plot$Freq.y * 100, 2)
  all.tissues.promoter.table.category.common.plot$Var1 = factor(all.tissues.promoter.table.category.common.plot$Var1, category.order)
  p2=ggplot(all.tissues.promoter.table.category.common.plot,aes(x=Var1,y=percentage,fill=Var1)) +
    geom_bar(stat="identity", alpha=0.8,width=0.7) + #,position="fill"
    geom_text(aes(label=rawvalue),position = position_stack(reverse = FALSE, vjust=0.5), size=3.5) +
    scale_fill_manual(values=c(rep("#ffac41", 2), rep("#77d8d8", 3), rep("#f76a8c", 3)) ) + guides(fill=FALSE) +
    xlab("") +
    ylab("Cateogry percentage (%)") +
    coord_flip() +
    ggtitle(paste(name,"common promoters")) +
    theme_classic.add + theme(plot.title = element_text(hjust = 0.5))
  print(p2)
  ggsave(p2,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.common.barplot.pdf")), width=3, height=3) #Fig. 4 left panel
  
  ##common in 2 tissues: 
  if(name %in% c("Brain","Cerebellum")){
    all.tissues.promoter.cluster.2common = c(rownames(all.tissues.promoter.up.cluster.2common), rownames(all.tissues.promoter.down.cluster.2common))
    all.tissues.promoter.cluster.2common.geneid = unique(all.tissues.promoter.cluster.category.selected[rownames(all.tissues.promoter.cluster.category.selected) %in% all.tissues.promoter.cluster.2common, "geneid"])
    all.tissues.promoter.table.category.2common = all.tissues.promoter.table.category.selected[rownames(all.tissues.promoter.table.category.selected) %in% all.tissues.promoter.cluster.2common.geneid, ]
    all.tissues.promoter.table.category.2common.plot = data.frame(table(all.tissues.promoter.table.category.2common$category))
    save(all.tissues.promoter.table.category.2common, file=file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.2common.RData")))
    
    all.tissues.promoter.table.category.2common.plot = merge(data.frame(table(all.tissues.promoter.table.category.2common$category)),table.category.statistics, by.x="Var1", by.y="Var1")
    all.tissues.promoter.table.category.2common.plot$rawvalue = all.tissues.promoter.table.category.2common.plot$Freq.x
    all.tissues.promoter.table.category.2common.plot$percentage = round(all.tissues.promoter.table.category.2common.plot$Freq.x/all.tissues.promoter.table.category.2common.plot$Freq.y * 100, 2)
    all.tissues.promoter.table.category.2common.plot$Var1 = factor(all.tissues.promoter.table.category.2common.plot$Var1, category.order)
    p2=ggplot(all.tissues.promoter.table.category.2common.plot,aes(x=Var1,y=percentage,fill=Var1)) +
      geom_bar(stat="identity", alpha=0.8,width=0.7) + #,position="fill"
      geom_text(aes(label=rawvalue),position = position_stack(reverse = FALSE, vjust=0.5), size=3.5) +
      scale_fill_manual(values=c(rep("#ffac41", 2), rep("#77d8d8", 3), rep("#f76a8c", 3)) ) + guides(fill=FALSE) +
      xlab("") +
      ylab("Cateogry percentage (%)") +
      coord_flip() +
      ggtitle(paste("Brain Cerebellum common promoters")) +
      theme_classic.add + theme(plot.title = element_text(hjust = 0.5))
    print(p2)
    ggsave(p2,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.2common.barplot.pdf")), width=3, height=3)
  }
  
  ##specific to one tissue
  all.tissues.promoter.up.cluster.specific = all.tissues.promoter.up.cluster.overlap[all.tissues.promoter.up.cluster.overlap[,name]==1, ]   # & rowSums(all.tissues.promoter.up.cluster.overlap) %in% 1
  all.tissues.promoter.down.cluster.specific = all.tissues.promoter.down.cluster.overlap[all.tissues.promoter.down.cluster.overlap[,name]==1, ]  # & rowSums(all.tissues.promoter.down.cluster.overlap) %in% 1
  if(name == "Brain" | name == "Cerebellum"){all.tissues.promoter.up.cluster.specific = rbind(all.tissues.promoter.up.cluster.specific, all.tissues.promoter.up.cluster.2common)}
  if(name == "Brain" | name == "Cerebellum"){all.tissues.promoter.down.cluster.specific = rbind(all.tissues.promoter.down.cluster.specific, all.tissues.promoter.down.cluster.2common)}
  all.tissues.promoter.cluster.specific = c(rownames(all.tissues.promoter.up.cluster.specific), rownames(all.tissues.promoter.down.cluster.specific))
  all.tissues.promoter.cluster.specific.geneid = unique(all.tissues.promoter.cluster.category.selected[rownames(all.tissues.promoter.cluster.category.selected) %in% all.tissues.promoter.cluster.specific, "geneid"])
  all.tissues.promoter.cluster.specific.geneid = all.tissues.promoter.cluster.specific.geneid[!all.tissues.promoter.cluster.specific.geneid %in% all.tissues.promoter.cluster.common.geneid]
  #if(name %in% c("Brain","Cerebellum")){all.tissues.promoter.cluster.specific.geneid = all.tissues.promoter.cluster.specific.geneid[!all.tissues.promoter.cluster.specific.geneid %in% all.tissues.promoter.cluster.2common.geneid]}
  all.tissues.promoter.table.category.specific = all.tissues.promoter.table.category.selected[rownames(all.tissues.promoter.table.category.selected) %in% all.tissues.promoter.cluster.specific.geneid, ]
  save(all.tissues.promoter.table.category.specific, file=file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.specific.RData")))
  
  all.tissues.promoter.table.category.specific.plot = merge(data.frame(table(all.tissues.promoter.table.category.specific$category)),table.category.statistics, by.x="Var1", by.y="Var1")
  all.tissues.promoter.table.category.specific.plot$rawvalue = all.tissues.promoter.table.category.specific.plot$Freq.x
  all.tissues.promoter.table.category.specific.plot$percentage = round(all.tissues.promoter.table.category.specific.plot$Freq.x/all.tissues.promoter.table.category.specific.plot$Freq.y * 100, 2)
  all.tissues.promoter.table.category.specific.plot$Var1 = factor(all.tissues.promoter.table.category.specific.plot$Var1, category.order)
  p1=ggplot(all.tissues.promoter.table.category.specific.plot,aes(x=Var1,y=percentage,fill=Var1)) +
    geom_bar(stat="identity", alpha=0.8,width=0.7) + #,position="fill"
    geom_text(aes(label=rawvalue),position = position_stack(reverse = FALSE, vjust=0.5), size=3.5) +
    scale_fill_manual(values=c(rep("#ffac41", 2), rep("#77d8d8", 3), rep("#f76a8c", 3)) ) + guides(fill=FALSE) +
    xlab("") +
    ylab("Cateogry percentage (%)") +
    coord_flip() +
    ggtitle(paste(name,"specific promoters")) +
    theme_classic.add + theme(plot.title = element_text(hjust = 0.5))
  print(p1)
  ggsave(p1,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.specific.barplot.pdf")), width=3, height=3) #Fig. 4 middle panel
  
  
  # all.tissues.promoter.up.cluster.overlap.brain.cere = all.tissues.promoter.up.cluster.overlap[all.tissues.promoter.up.cluster.overlap$Brain==1 & all.tissues.promoter.up.cluster.overlap$Cerebellum==1 & rowSums(all.tissues.promoter.up.cluster.overlap[,-1:-2]) == 0, ]
  # all.tissues.promoter.down.cluster.overlap.brain.cere = all.tissues.promoter.down.cluster.overlap[all.tissues.promoter.down.cluster.overlap$Brain==1 & all.tissues.promoter.down.cluster.overlap$Cerebellum==1 & rowSums(all.tissues.promoter.down.cluster.overlap[,-1:-2]) == 0, ]
  # all.tissues.promoter.up.cluster.overlap.common = all.tissues.promoter.up.cluster.overlap[rowSums(all.tissues.promoter.up.cluster.overlap) >=2 & !(rownames(all.tissues.promoter.up.cluster.overlap) %in% rownames(all.tissues.promoter.up.cluster.overlap.brain.cere)) , ]
  # all.tissues.promoter.down.cluster.overlap.common = all.tissues.promoter.down.cluster.overlap[rowSums(all.tissues.promoter.down.cluster.overlap) >=2 & 
  #                                                                                             !(rownames(all.tissues.promoter.down.cluster.overlap) %in% rownames(all.tissues.promoter.down.cluster.overlap.brain.cere)) &
  #                                                                                               !(rownames(all.tissues.promoter.down.cluster.overlap) %in% rownames(all.tissues.promoter.down.cluster.overlap.testis.liver)), ]
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}



################################################
##GO annotation of overlapped genes and tissue specific genes
################################################
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(DOSE)

for (name in tissue.order){
  if(name=="Ovary"){gender="Female"} else {gender="Male"}
  # name = "Testis" #"Brain" 
  Pcutoff=0.05
  Qcutoff=0.05
  
  #load common and specific DE promoter list
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"gene.expression.star.mean.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.common.RData")))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.specific.RData")))
  if(name %in% c("Brain","Cerebellum")){load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.2common.RData")))}
  ###covert GeneID to entrez ID
  database.name = ifelse(length(grep("Mouse",prewd))!=0, "org.Mm.eg.db", "org.Hs.eg.db")
  species.name = ifelse(length(grep("Mouse",prewd))!=0, "Mus musculus", "Homo sapiens")
  expressed.genes.entrezid = bitr(row.names(gene.expression.star.mean[rowSums(gene.expression.star.mean)>10,]), fromType="ENSEMBL", toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=database.name)
  all.genes.entrezid = bitr(row.names(gene.expression.star.mean), fromType="ENSEMBL", toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=database.name)
  
  ##convert common and specific DE promoter list ###
  DEgenes.name.list = list(common = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.common),"ENTREZID"], 
                           specific = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.specific),"ENTREZID"])
  # if(name %in% c("Brain","Cerebellum")){
  #   DEgenes.name.list = list(common = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.common),"ENTREZID"],
  #                            specific = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.specific),"ENTREZID"],
  #                            common2 = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.2common),"ENTREZID"] )
  # }
  #####################
  ##GO Analysis
  #####################
  if(!file.exists(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.go.RData")))){ #TRUE
    all.tissues.promoter.table.category.go = try(compareCluster(DEgenes.name.list, fun = "enrichGO", OrgDb = database.name, ont="BP", pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID, readable = TRUE),silent = TRUE)
    save(all.tissues.promoter.table.category.go, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.go.RData"))))
  }
  #####dotPlot of enrich annotations
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.go.RData")))
  pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.go.pdf") ), width = 5, height = 3) #Fig. 4 right panel
  print(dotplot(all.tissues.promoter.table.category.go, showCategory=5))
  dev.off()
  
  #####################
  ##MSigDb C5 (Not used anymore)
  #####################
  if(FALSE){
    #####common DE promoter
    all.tissues.promoter.table.category.common.msigc5 = try(enricher(DEgenes.name.list[[1]], TERM2GENE=msigdbr(species = species.name , category = "C5") %>% dplyr::select(gs_name, entrez_gene), pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID),silent = TRUE)
    save(all.tissues.promoter.table.category.common.msigc5, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.common.msigc5.RData"))))
    #####specific DE promoter
    all.tissues.promoter.table.category.specific.msigc5 = try(enricher(DEgenes.name.list[[2]], TERM2GENE=msigdbr(species = species.name , category = "C5") %>% dplyr::select(gs_name, entrez_gene), pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID),silent = TRUE)
    save(all.tissues.promoter.table.category.specific.msigc5, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.specific.msigc5.RData"))))
    #####DE promoter that are common in "Brain","Cerebellum" 
    if(name %in% c("Brain","Cerebellum")){
      all.tissues.promoter.table.category.2common.msigc5 = try(enricher(DEgenes.name.list[[3]], TERM2GENE=msigdbr(species = species.name , category = "C5") %>% dplyr::select(gs_name, entrez_gene), pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID),silent = TRUE)
      save(all.tissues.promoter.table.category.2common.msigc5, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.2common.msigc5.RData"))))
    }
    #####dotPlot of enrich annotations
    pdf( file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.msigc5.pdf") ))
    print(dotplot(all.tissues.promoter.table.category.common.msigc5, showCategory=14))
    if(name %in% c("Brain","Cerebellum")){print(dotplot(all.tissues.promoter.table.category.2common.msigc5, showCategory=14))}
    print(dotplot(all.tissues.promoter.table.category.specific.msigc5, showCategory=14))
    dev.off()
  }
  
  if(length(grep("Mouse",prewd))==0){
    #####################
    ### Disease Analysis ###
    #####################
    ##DO
    all.tissues.promoter.table.category.do = try(compareCluster(DEgenes.name.list, fun = "enrichDO", ont="DO", pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID, readable = TRUE),silent = TRUE)
    save(all.tissues.promoter.table.category.do, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.do.RData"))))
    #####dotPlot of enrich annotations
    pdf( file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.do.pdf") ))
    print(dotplot(all.tissues.promoter.table.category.do, showCategory=14))
    dev.off()
  }
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

###GO annotation comparision between all tissues  (Revision)
if(TRUE){
  all.tissues.promoter.table.category.go.list <- list()
  for (name in tissue.order){
    if(name=="Ovary"){gender="Female"} else {gender="Male"}
    df = as.data.frame(get(load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.go.RData")))))
    df$Tissue = name
    if (name %in% c("Heart")){df=df[df$qvalue<0.01,]}
    if (name %in% c("Brain")){df=as.data.frame(rbind(df[df$qvalue<0.05 & df$Cluster=="common",],df[df$qvalue<0.01 & df$Cluster=="specific",]))}
    all.tissues.promoter.table.category.go.list[[name]] = df
  }
  all.tissues.promoter.table.category.go = do.call(rbind, all.tissues.promoter.table.category.go.list)  
  save(all.tissues.promoter.table.category.go, file=file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",paste0("all.tissues.promoter.table.category.go.upsetr.RData") ))
  write.table(all.tissues.promoter.table.category.go, file=file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",paste0("all.tissues.promoter.table.category.go.upsetr.txt") ), quote = FALSE, row.names = FALSE, col.names = TRUE) #Supplementary table S4 in the paper
  
  # Load required library
  build.plot.upsetr = function(upset_data,filename){
    library(UpSetR)
    pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",paste0("all.tissues.promoter.table.category.go.upsetr.", filename,".pdf") ), width = 5, height = 4.5, onefile = FALSE) #Fig. S9a-b
    print(upset(upset_data, 
          order.by = "freq",         # Order intersections by frequency
          main.bar.color = "steelblue", 
          sets.bar.color = "darkorange",
          keep.order = TRUE, 
          text.scale = 1.5,
          nsets = 7,                  # Number of sets to display
          nintersects = NA) )         # Scale text for better readability
    dev.off()
  }
  
  # For common GO terms
  all.tissues.promoter.table.category.go.qvalue.common <- subset(all.tissues.promoter.table.category.go, Cluster == "common")
  all.tissues.promoter.table.category.go.qvalue.common.selected <- split(all.tissues.promoter.table.category.go.qvalue.common$ID, all.tissues.promoter.table.category.go.qvalue.common$Tissue)
  all.tissues.promoter.table.category.go.qvalue.common.upset <- fromList(all.tissues.promoter.table.category.go.qvalue.common.selected)
  build.plot.upsetr(all.tissues.promoter.table.category.go.qvalue.common.upset, filename="common")
  
  # For specific GO terms
  all.tissues.promoter.table.category.go.qvalue.specific <- subset(all.tissues.promoter.table.category.go, qvalue < 0.05 & Cluster == "specific")
  all.tissues.promoter.table.category.go.qvalue.specific.selected <- split(all.tissues.promoter.table.category.go.qvalue.specific$ID, all.tissues.promoter.table.category.go.qvalue.specific$Tissue)
  all.tissues.promoter.table.category.go.qvalue.specific.upset <- fromList(all.tissues.promoter.table.category.go.qvalue.specific.selected)
  build.plot.upsetr(all.tissues.promoter.table.category.go.qvalue.specific.upset, filename="specific")
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

################################################
##TF Motif analysis of overlapped genes and tissue specific genes (Revision)
################################################
if (TRUE){
  library(BSgenome.Mmusculus.UCSC.mm10) # Mouse genome
  library(TFBSTools)
  library(JASPAR2020)
  library(motifmatchr)
  library(PWMEnrich)
  library("PWMEnrich.Mmusculus.background")
  data("PWMLogn.mm9.MotifDb.Mmus")
  
  ##promoter_annotation
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"absolutePromoterActivityTable.RData"))
  absolutePromoterActivityTable = absolutePromoterActivityTable[,1:7]
  absolutePromoterActivityTableGranges = makeGRangesFromDataFrame(absolutePromoterActivityTable, keep.extra.columns = TRUE)
  
  ##significant time changing isoforms in all tissues
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration", "all.tissues.promoter.table.category.list.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration", "all.tissues.promoter.cluster.category.list.RData"))
  
  for (name in tissue.order){ #"Brain"
    if(name=="Ovary"){gender="Female"} else {gender="Male"}
    # Step 1: Extract the promoter data for major and minor promoters
    promoter_data <- all.tissues.promoter.table.category.list[[name]]
    tf.motif.enrichment = rbind()
    tf.motif.enrichment.list = list()
    
    for (promoter.category in unique(promoter_data$category)){
      #promoter.category="Up_Up"
      print(promoter.category)
      # Separate Major and Minor promoters
      major_promoters <- as.character(promoter_data[promoter_data$category == promoter.category,"Major"])
      minor_promoters <- as.character(promoter_data[promoter_data$category == promoter.category,"Minor"])
      minor_promoters <- unlist(strsplit(minor_promoters, split = "_"))
      
      
      # Step 2: Map promoters to coordinates in absolutePromoterActivityTableGranges
      # Extract major and minor promoter coordinates
      major_coords <- promoters(absolutePromoterActivityTableGranges[absolutePromoterActivityTableGranges$promoterId %in% major_promoters], upstream = 200, downstream = 0)
      minor_coords <- promoters(absolutePromoterActivityTableGranges[absolutePromoterActivityTableGranges$promoterId %in% minor_promoters], upstream = 200, downstream = 0)
      
      # Step 3: Get sequences for promoter regions
      major_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, major_coords)
      minor_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, minor_coords)
      
      # Step 4: Compute motif enrichment for major and minor sequences
      # Major sequences
      major_motif_scores <- motifEnrichment(major_seqs, PWMLogn.mm9.MotifDb.Mmus)
      major_motif_results = groupReport(major_motif_scores)
      
      # Minor sequences
      minor_motif_scores <- motifEnrichment(minor_seqs, PWMLogn.mm9.MotifDb.Mmus)
      minor_motif_results = groupReport(minor_motif_scores)
      
      # Step 5: save all the enriched motifs in a dataframe
      tf.motif.enrichment <- rbind(tf.motif.enrichment,
        cbind(as.data.frame(major_motif_results[1:100,]), group = "Major", promoter_category = promoter.category),
        cbind(as.data.frame(minor_motif_results[1:100,]), group = "Minor", promoter_category = promoter.category)
      )
      tf.motif.enrichment.list[[promoter.category]] = list(major=major_motif_results, minor=minor_motif_results)
            
    }
    
    # Save results to file
    save(tf.motif.enrichment.list, file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.motif.RData")) )
    write.table(tf.motif.enrichment, file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.motif.txt")), sep = "\t", row.names = FALSE, quote = FALSE)  
  }
  #cat Mammalian2019Promoter/CM2019RNAseqMouse/Output/StageSpecific/Proactiv/edger/*/*/all.tissues.promoter.table.category.motif.txt > Mammalian2019Promoter/CM2019RNAseqMouse/Output/StageSpecific/Proactiv/edger/Integration/all.tissues.promoter.table.category.motif.txt #Supplementary table S4 in the paper
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

##Plot
for (name in tissue.order){
  print(name)
  if(name=="Ovary"){gender="Female"} else {gender="Male"}
  load(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.motif.RData")))
  
  for (promoter.category in names(tf.motif.enrichment.list)){
    #promoter.category="Up_Up"
    # print(promoter.category )
    major_motif_results = tf.motif.enrichment.list[[promoter.category]][[1]]
    minor_motif_results = tf.motif.enrichment.list[[promoter.category]][[2]]
    
    ##only extract p.value < 0.05
    if(promoter.category %in% c("Down_Up", "Up_Down")){cutoff.number = 1e-2} else {cutoff.number = 1e-2}
    major_motif_results = major_motif_results[major_motif_results$p.value < cutoff.number & major_motif_results$top.motif.prop > 0.1,]
    minor_motif_results = minor_motif_results[minor_motif_results$p.value < cutoff.number & minor_motif_results$top.motif.prop > 0.1,]
    
    # ##order based on other columns
    # major_motif_results = major_motif_results[order(major_motif_results$top.motif.prop, decreasing = T),]
    # minor_motif_results = minor_motif_results[order(minor_motif_results$top.motif.prop, decreasing = T),]
    
    ##Top 10 in either Major or Minor sequences, top.motif.prop means the proportion of sequences owns this motif
    top10.sequences = unique(c(major_motif_results$id, rev(minor_motif_results$id)))
    # print(paste0("Overlap between Major and Minor:", overlap <- length(intersect(major_motif_results$id, minor_motif_results$id))/length(top10.sequences)))
    major_motif_results_top10 = as.data.frame(major_motif_results)[major_motif_results$id %in% top10.sequences, ]
    minor_motif_results_top10 = as.data.frame(minor_motif_results)[minor_motif_results$id %in% top10.sequences, ]
    
    # Step 5: Prepare data for dot plot visualization
    major_motif_results_top10$group="Major"
    minor_motif_results_top10$group="Minor"
    plot_df <- rbind(major_motif_results_top10, minor_motif_results_top10)
    plot_df$target = factor(plot_df$target, levels = rev(unique(plot_df$target)))
    
    # Step 6: Plot the results (dotplot)
    p = ggplot(plot_df, aes(x = group, y = target, color = top.motif.prop*100, size = -log10(p.value))) + # raw.score
      geom_point() +  # Use fill for an additional variable aes(fill = top.motif.prop), shape = 21
      scale_color_gradient2(
        low = "#FFF078",
        mid = "#D91656",
        high = "#4F1787",
        midpoint = median(-log10(plot_df$p.value), na.rm = TRUE)  ) +
      # scale_color_gradient(
      #   low = "#D91656", 
      #   high = "#4F1787") + 
      theme_minimal() +
      labs(
        title = "",
        x = "",
        y = "",
        size = "-log10(P value)",
        color = "Motif frequency"
        # fill = "Motif Proportion"
      ) +
      theme(axis.text.y = element_text(angle = 0, hjust = 1),axis.text.x = element_text(angle = 30), axis.text = element_text(size = 14, colour = "black", family = "Helvetica"), legend.position = "bottom",legend.box = "horizontal") +  
      # guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
      scale_size_continuous(range = c(3, 8))  # Adjust dot sizes
    # print(p)
    ggsave(plot=p, filename=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.motifdot.", promoter.category,".pdf") ), width = 2, height = 10) #Fig. 4
    
    # Step 7: Plot the results (motif sequence)
    pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.motifseq.", promoter.category,".pdf")), width = 5, height = 5)
    par(mfrow = c(2, 1))
    plot(major_motif_results, fontsize = 7, id.fontsize = 5)
    plot(minor_motif_results, fontsize = 7, id.fontsize = 5)
    dev.off()
  }
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}
  

