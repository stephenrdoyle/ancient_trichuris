## ANGSD
- exploring the use of ANGSD, whcih can use genotype likelihoods for a number of analyses. Probably good for the low coverage datasets
- started off calculating "identity by state" (IBS), which is another way of showing genetic similarity between samples

```bash

cd /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/05_ANALYSIS/ANGSD

bsub.py --threads 5 10 angsd "/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/angsd/angsd -bam bam.list -minMapQ 30 -minQ 20 -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -minMaf 0.05 -P 5"

```

```R
# load libraries
library(gplots)

# load data
cov_data <- as.matrix(read.table("angsdput.covMat",header=F))
ibs_data <- as.matrix(read.table("angsdput.ibsMat",header=F))

names <- as.matrix(read.table("sample.names"))
rownames(ibs_data) <- names
colnames(ibs_data) <- names
rownames(cov_data) <- names
colnames(cov_data) <- names


pdf("nuclear_ibs_heatmap.pdf",width=15, height=15)
heatmap.2(ibs_data, trace="none", margins=c(12,12))
dev.off()

pdf("nuclear_covariance_heatmap.pdf",width=15, height=15)
cov_data[cov_data > 1] <- 1
heatmap.2(cov_data, trace="none", margins=c(12,12))
dev.off()

```

Figure: [IBS for nuclear markers](../04_analysis/nuclear_ibs_heatmap.pdf)

Figure: [covariance of IBS for nuclear markers](../04_analysis/nuclear_covariance_heatmap.pdf)
