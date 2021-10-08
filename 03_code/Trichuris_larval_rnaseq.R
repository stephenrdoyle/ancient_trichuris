# Analysis of Trichuris larvae bulk RNA-seq

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(ggplot2)
options(stringsAsFactors = FALSE)
setwd('/Users/ar11/Work/Writing/Trichuris/Larval RNA-seq/Trichuris_Larval_RNAseq_R')
workdir = '/Users/ar11/Work/Writing/Trichuris/Larval RNA-seq/Trichuris_Larval_RNAseq_R/'

countData<- read.table("gene_counts.dat", sep = "\t", header = TRUE, row.names=1)
# Remove bleached samples
countData <- countData[,colnames(countData) != c("Bleach_Hatched_Egg7.10", "Bleach_Hatched_Egg8.20", "Bleach_Hatched_Egg9.30")]

colData<-read.table("samples_full.dat", sep = "\t", header = TRUE, row.names=1)
# Remove bleached samples
colData <- colData[rownames(colData) != c("Bleach_Hatched_Egg7.10", "Bleach_Hatched_Egg8.20", "Bleach_Hatched_Egg9.30"),]

# Reorder count columns to match colData rows
countData <- countData[rownames(colData)]

desc <- read.table("trichuris_muris.PRJEB126.WBPS14.gene.desc.txt", header=TRUE, row.names=1, sep="\t")
pfam <- read.table("trichuris_muris.PRJEB126.WBPS14.gene.pfam.txt", header=TRUE, sep="\t")
goterms <- read.table("trichuris_muris.PRJEB126.WBPS14.gene.go.txt", header=FALSE, sep="\t")
signalp_tmhmm_data <- read.table("signalp_tmhmm_wormbaseparasite.txt", header=TRUE, sep='\t')

# Process SignalP/TMHMM data
#signalp = subset(signalp_tmhmm_data, signalp_tmhmm_data$SignalP. == "SignalP-noTM")$Gene.name
#tmhmm = subset(signalp_tmhmm_data, signalp_tmhmm_data$TM == "TMhelix")$Gene.name
rownames(signalp_tmhmm_data) = signalp_tmhmm_data$Gene.name
signalp_tmhmm_data$Gene.name = c()
signalp_tmhmm_data$TM = c()
#########################
# DE and GO enrichment
#
# n.b. python script used to convert 
# DE results into Excel workbook
#########################
dds <- DESeqDataSetFromMatrix(countData= countData,colData= colData, design=~condition)
# Filter samples/genes
dds<-  dds[    rowSums( counts(dds)) > 1, ]
dds<-  dds[    ,colSums( counts(dds)) > 500000 ] # Filter samples with low read counts
dds<-DESeq(dds)

# Compare Egg and L1_3h
res_egg_L1_3hr <- results(dds, contrast=c("condition", "L1_3hr", "Unhatched"))
 #summary(res_egg_L1_3hr)
sum(res_egg_L1_3hr$padj < 0.01, na.rm=TRUE)
down_egg_L1_3hr <- rownames(subset(res_egg_L1_3hr, padj < 0.05 & log2FoldChange < -1))
up_egg_L1_3hr <- rownames(subset(res_egg_L1_3hr, padj < 0.05 & log2FoldChange > 1))
runtopGO(down_egg_L1_3hr, "trichuris_muris.PRJEB126.WBPS14.gene.go.txt", "down_egg_L1_3hr", workdir, FDR=0.05)
runtopGO(up_egg_L1_3hr, "trichuris_muris.PRJEB126.WBPS14.gene.go.txt", "up_egg_L1_3hr", workdir, FDR=0.05)

# Compare L1_3h and L1_24h
res_L1_3h_L1_24hr <- results(dds, contrast=c("condition","L1_24hr","L1_3hr"))
summary(res_L1_3h_L1_24hr)
sum(res_L1_3h_L1_24hr$padj < 0.01, na.rm=TRUE)
down_L1_3h_L1_24hr <- rownames(subset(res_L1_3h_L1_24hr, padj < 0.05 & log2FoldChange < -1))
up_L1_3h_L1_24hr <- rownames(subset(res_L1_3h_L1_24hr, padj < 0.05 & log2FoldChange > 1))
runtopGO(down_L1_3h_L1_24hr, "trichuris_muris.PRJEB126.WBPS14.gene.go.txt", "down_L1_3h_L1_24hr", workdir, FDR=0.05)
runtopGO(up_L1_3h_L1_24hr, "trichuris_muris.PRJEB126.WBPS14.gene.go.txt", "up_L1_3h_L1_24hr", workdir, FDR=0.05)

# Compare L1_24h and L2
res_L1_24hr_L2 <- results(dds, contrast=c("condition","L2","L1_24hr"))
summary(res_L1_24hr_L2)
sum(res_L1_24hr_L2$padj < 0.01, na.rm=TRUE)
down_L1_24hr_L2 <- rownames(subset(res_L1_24hr_L2, padj < 0.05 & log2FoldChange < -1))
up_L1_24hr_L2 <- rownames(subset(res_L1_24hr_L2, padj < 0.05 & log2FoldChange > 1))
runtopGO(down_L1_24hr_L2, "trichuris_muris.PRJEB126.WBPS14.gene.go.txt", "down_L1_24hr_L2", workdir, FDR=0.05)
runtopGO(up_L1_24hr_L2, "trichuris_muris.PRJEB126.WBPS14.gene.go.txt", "up_L1_24hr_L2", workdir, FDR=0.05)

#####################
# Draw/Write out PCA
#####################
rld<- rlog(dds,blind=  FALSE)
condition_pca <- make_pca(rld, "condition"); condition_pca

cairo_pdf("early_larvae_pca.pdf", width=6, height=3)
condition_pca
dev.off()

################
# FPKM heatmaps
################
##### combine read counts for replicates
egg_counts <- rowSums(countData[,c("Unhatched_LB_Egg7.8", "Unhatched_Egg8.17", "Unhatched_Egg9.27",
             "Unhatched_LB_Egg8.18", "Unhatched_LB_Egg9.28")])
l1_3hr_counts <- rowSums(countData[,c("L1_3hr_mus1_Egg8.13", "L1_3hr_mus1_Egg9.23", "L1_3hr_mus2_Egg7.4")])
l1_24hr_counts <-rowSums(countData[,c("L1_24hr_mus2_Egg7.6", "L1_24hr_mus1_Egg8.15", "L1_24hr_mus1_Egg9.25",
                     "L1_24hr_mus2_Egg8.16", "L1_24hr_mus2_Egg9.26")])
l2_counts <- rowSums(countData[,c("L2_mus1_Egg7.1", "L2_mus2_Egg7.2", "L2_mus1_Egg8.11", "L2_mus1_Egg9.21",
                                  "L2_mus2_Egg8.12", "L2_mus2_Egg9.22")])

combined_counts <- data.frame(row.names=rownames(countData), egg=egg_counts, l1_3hr=l1_3hr_counts,
                              l1_24hr=l1_24hr_counts, l2=l2_counts)

##### Calculate FPKMs
 # gene lengths (length of longest transcript for each gene)
lengths <- read.table("trichuris_muris.PRJEB126.WBPS14.mRNA_transcripts.lengths.longest", header=FALSE, row.names=1)
 # library sizes
egg_lib <- sum(egg_counts)
l1_3hr_lib <- sum(l1_3hr_counts)
l1_24hr_lib <- sum(l1_24hr_counts)
l2_lib <- sum(l2_counts)

combined_fpkm = combined_counts # copy counts dataframe to be overwriten with FPKMs

##### Loop over genes in counts dataframe to make FPKM dataframe
for (i in 1:nrow(combined_counts)){
  # calculate FPKM for each colmn
  egg_fpkm = (combined_counts[i,]$egg * 10^9) / (egg_lib * lengths[rownames(combined_counts[i,]),])
  l1_3hr_fpkm = (combined_counts[i,]$l1_3hr * 10^9) / (l1_3hr_lib * lengths[rownames(combined_counts[i,]),])
  l1_24hr_fpkm = (combined_counts[i,]$l1_24hr * 10^9) / (l1_24hr_lib * lengths[rownames(combined_counts[i,]),])
  l2_fpkm = (combined_counts[i,]$l2 * 10^9) / (l2_lib * lengths[rownames(combined_counts[i,]),])
 
  # add FPKMs for each column to fpkm dataframe
  combined_fpkm[i,]$egg = egg_fpkm
  combined_fpkm[i,]$l1_3hr = l1_3hr_fpkm
  combined_fpkm[i,]$l1_24hr = l1_24hr_fpkm
  combined_fpkm[i,]$l2 = l2_fpkm
}

# get list of all DE genes
all_de <- unique(c(down_egg_L1_3hr, up_egg_L1_3hr, down_L1_3h_L1_24hr, up_L1_3h_L1_24hr, down_L1_24hr_L2, up_L1_24hr_L2))

##### WAPs heatmap
# identify WAP domain-containing genes by Pfam domain
waps = subset(pfam, pfam == "PF00095")$Gene_name
# extract expression patterns for WAPs which are DE
waps_fpkm = subset(combined_fpkm, rownames(combined_fpkm) %in% waps & rownames(combined_fpkm) %in% all_de)
# Draw heatmap
pheatmap(log(waps_fpkm + 1, base=2), cluster_cols=FALSE, labels_col = c("Egg", "L1 3h", "L1 24h", "L2"),
         annotation_row = signalp_tmhmm_data)
pheatmap(log(waps_fpkm + 1, base=2), cluster_cols=FALSE, labels_col = c("Egg", "L1 3h", "L1 24h", "L2"),
         annotation_row = signalp_tmhmm_data, filename="waps_heatmap.pdf")

##### Serine endopeptidases heatmap
# Identify serine type endopeptidase by GO term
ste <- subset(goterms, grepl("GO:0004252", goterms$V2))$V1
# extract expression patterns for serine endopeptidases
ste_fpkm = subset(combined_fpkm, rownames(combined_fpkm) %in% ste & rownames(combined_fpkm) %in% all_de)
# Draw heatmap
pheatmap(log(ste_fpkm + 1, base=2), cluster_cols=FALSE, labels_col = c("Egg", "L1 3h", "L1 24h", "L2"),
         annotation_row = signalp_tmhmm_data)
pheatmap(log(ste_fpkm + 1, base=2), cluster_cols=FALSE, labels_col = c("Egg", "L1 3h", "L1 24h", "L2"),
         annotation_row = signalp_tmhmm_data, filename="ste_heatmap.pdf")

dev.off()

##### Make a file of Peptidase inhibitors which are DE
pepinhib <- subset(goterms, grepl("GO:0030414", goterms$V2))$V1
pepinhib_fpkm = subset(combined_fpkm, rownames(combined_fpkm) %in% pepinhib & rownames(combined_fpkm) %in% all_de)
write.table(rownames(pepinhib_fpkm), file="pepinhib_de.list", quote=FALSE, sep="\t", 
            row.names=FALSE, col.names=FALSE)

##### Kunitz domain protease inhibitors heatmap
kunitz <- subset(pfam, pfam == "PF00014")$Gene_name
# extract expression patterns for WAPs which are DE
kunitz_fpkm = subset(combined_fpkm, rownames(combined_fpkm) %in% kunitz & rownames(combined_fpkm) %in% all_de)
# Draw heatmap
pheatmap(log(kunitz_fpkm + 1, base=2), cluster_cols=FALSE, labels_col = c("Egg", "L1 3h", "L1 24h", "L2"),
         annotation_row = signalp_tmhmm_data)
pheatmap(log(kunitz_fpkm + 1, base=2), cluster_cols=FALSE, labels_col = c("Egg", "L1 3h", "L1 24h", "L2"),
         annotation_row = signalp_tmhmm_data, filename="kunitz_heatmap.pdf")

###################################
##### FPKMs for separate replicates
###################################
separate_fpkm = countData # copy counts dataframe to be overwriten with FPKMs

##### Loop over genes in counts dataframe to make FPKM dataframe
for (i in 1:nrow(countData)){
  #print(i)
  for (j in 1:ncol(countData)) {
    #print(j)
    fpkm_val = (countData[i,j] * 10^9) / (colSums(countData[j]) * lengths[rownames(countData[i,]),])
    separate_fpkm[i,j] = fpkm_val
  }
}

######
#separate heatmaps
# extract expression patterns for serine endopeptidases
ste_fpkm_separate = subset(separate_fpkm[c("Unhatched_LB_Egg7.8", "Unhatched_Egg8.17", "Unhatched_Egg9.27",
                                           "Unhatched_LB_Egg8.18", "Unhatched_LB_Egg9.28",
                                           "L1_3hr_mus1_Egg8.13", "L1_3hr_mus1_Egg9.23", "L1_3hr_mus2_Egg7.4",
                                           "L1_24hr_mus2_Egg7.6", "L1_24hr_mus1_Egg8.15", "L1_24hr_mus1_Egg9.25",
                                           "L1_24hr_mus2_Egg8.16", "L1_24hr_mus2_Egg9.26",
                                           "L2_mus1_Egg7.1", "L2_mus2_Egg7.2", "L2_mus1_Egg8.11", "L2_mus1_Egg9.21",
                                           "L2_mus2_Egg8.12", "L2_mus2_Egg9.22" )], rownames(separate_fpkm) %in% ste & rownames(separate_fpkm) %in% all_de)
# Draw heatmap
pheatmap(log(ste_fpkm_separate + 1, base=2), cluster_cols=FALSE)


#################################
# rld heatmaps (DEPRECATED!!!!)
#################################
# Subset samples for relevant conditions
rld_heatmap = rld[,c("Unhatched_LB_Egg7.8", "Unhatched_Egg8.17", "Unhatched_Egg9.27",
                     "Unhatched_LB_Egg8.18", "Unhatched_LB_Egg9.28",
                     "L1_3hr_mus1_Egg8.13", "L1_3hr_mus1_Egg9.23", "L1_3hr_mus2_Egg7.4",
                     "L1_24hr_mus2_Egg7.6", "L1_24hr_mus1_Egg8.15", "L1_24hr_mus1_Egg9.25",
                     "L1_24hr_mus2_Egg8.16", "L1_24hr_mus2_Egg9.26",
                     "L2_mus1_Egg7.1", "L2_mus2_Egg7.2", "L2_mus1_Egg8.11", "L2_mus1_Egg9.21",
                     "L2_mus2_Egg8.12", "L2_mus2_Egg9.22" )]

# WAPs which are DE
waps_rld = subset(rld_heatmap, rownames(rld_heatmap) %in% waps & rownames(rld_heatmap) %in% all_de)
pheatmap(assay(waps_rld), cluster_cols=FALSE)

# Serine type endopeptidase by GO term
ste <- subset(goterms, grepl("GO:0004252", goterms$V2))$V1
ste_rld = subset(rld_heatmap, rownames(rld_heatmap) %in% ste & rownames(rld_heatmap) %in% all_de)
pheatmap(assay(ste_rld), cluster_cols=FALSE)

###################
# FUNCTIONS
###################
make_pca <- function (rld, design) {
  
  pcaData <- plotPCA(rld,intgroup=  c(design), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=get(design))) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  return (pca_plot)
  #return (pcaData)
}

##################
#TopGO function
###################

runtopGO <-function(x, go_file, outstem, outdir, FDR=0.01){
  library(topGO)  
  
  node_size <- 5
  method = 'weight01'
  
  go_file <- read.table(go_file, header=FALSE, stringsAsFactors=FALSE)
  names(go_file) = c('id', 'go')
  ref.vec <- strsplit(go_file$go, split=",", fixed=T)
  names(ref.vec) <- go_file$id
  all.ids <- go_file$id
  
  scores <- rep(0, nrow(go_file))
  names(scores) <- go_file$id
  scores[go_file$id %in% x] <- 1
  
  geneSelectionFun <- function(score){
    return(score >= 1)
  }
  
  GOdataBP <- new("topGOdata",  ontology = 'BP', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = node_size, description = '')
  GOdataMF <- new("topGOdata",  ontology = 'MF', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = node_size, description = '')
  GOdataCC <- new("topGOdata",  ontology = 'CC', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = node_size, description = '')
  
  resultTopgoBP <- runTest(GOdataBP,algorithm=method,statistic="Fisher")
  resultTopgoMF <- runTest(GOdataMF,algorithm=method,statistic="Fisher")
  resultTopgoCC <- runTest(GOdataCC,algorithm=method,statistic="Fisher")
  
  resBP<-GenTable( GOdataBP, topGO = resultTopgoBP, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)
  resMF<-GenTable( GOdataMF, topGO = resultTopgoMF, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)
  resCC<-GenTable( GOdataCC, topGO = resultTopgoCC, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)
  
  #outfile <- paste("/Users/ar11/trichuris/single_cell_infection/10x/10x_reps/cluster", clustern,"_go_bp.txt", sep="")
  outfile <- paste(outdir, outstem, "_go_bp.txt", sep="")
  print(outfile)
  write.table(resBP[as.numeric(resBP$topGO) < FDR,], file=outfile, sep="\t", quote=FALSE)
  outfile <- paste(outdir, outstem, "_go_mf.txt", sep="")
  print(outfile)
  write.table(resMF[as.numeric(resMF$topGO) < FDR,], file=outfile, sep="\t", quote=FALSE)
  outfile <- paste(outdir, outstem, "_go_cc.txt", sep="")
  print(outfile)
  write.table(resCC[as.numeric(resCC$topGO) < FDR,], file=outfile, sep="\t", quote=FALSE)
  
}
