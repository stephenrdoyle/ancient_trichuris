# Other analyses performed, but didnt make it into the paper - keeping them for a rainy day


## F stats using treemix
- ended up using q3stats rather than this
- could run 4 pop analyses, but didnt make sense with so few populations

```bash
threepop -i treemix.LDpruned.treemix.frq.gz -k 500 > threepop.out
fourpop -i treemix.LDpruned.treemix.frq.gz -k 500 > fourpop.out

cat threepop.out | grep ";" | grep -v "COLOBUS" | grep -v "LEAF" > threepop.out_2
cat fourpop.out | grep ";" | grep -v "COLOBUS" | grep -v "LEAF" > fourpop.out_2
```
- plots
```R
# load libraries
library(tidyverse)

# load data
data <- read.table("threepop.out_2",header=F)
colnames(data) <- c("three_populations", "f_3", "std_error", "z_score")

# plot f_3 (sorted) for each three population test
#-
ggplot(data, aes(f_3, reorder(three_populations, -f_3), col = z_score)) +
     geom_point(size = 2) +
     geom_segment(aes(x = f_3-std_error, y = three_populations, xend = f_3+std_error, yend = three_populations)) +
     theme_bw() +
     labs(x = "f3 statistic" , y = "")

ggsave("plot_treemix_f3_statistics.png")


# load data
data <- read.table("fourpop.out_2",header=F)
colnames(data) <- c("four_populations", "f_4", "std_error", "z_score")

# plot f_4 (sorted) for each three population test
#-
ggplot(data, aes(f_4, reorder(four_populations, -f_4), col = z_score)) +
     geom_point(size = 2) +
     geom_segment(aes(x = f_4-std_error, y = four_populations, xend = f_4+std_error, yend = four_populations)) +
     theme_bw() +
     labs(x = "f4 statistic" , y = "")

ggsave("plot_treemix_f4_statistics.png")

```
- f3 stats
![](../04_analysis/plot_treemix_f3_statistics.png)

- f4 stats
![](../04_analysis/plot_treemix_f4_statistics.png)









### Check depth per genotype per sample, and mask low coverage gts
```bash
bsub.py 1 select_mitoSNPs_GT \
"gatk VariantsToTable \
  --reference ${REFERENCE} \
  --variant ${VCF%.vcf.gz}.mitoALL.filtered.vcf \
  --fields CHROM --fields POS --genotype-fields GT --genotype-fields DP \
  --output ${VCF%.vcf.gz}.mitoALL.filtered.DP.table"

 bsub.py 1 select_nuclearSNPs_GT \
 "gatk VariantsToTable \
  --reference ${REFERENCE} \
  --variant ${VCF%.vcf.gz}.nuclearALL.filtered.vcf \
  --fields CHROM --fields POS --genotype-fields GT --genotype-fields DP \
  --output ${VCF%.vcf.gz}.nuclearALL.filtered.DP.table"



# make per sample depth datasets
ncol=$(awk 'NR==1{print NF}' ${VCF%.vcf.gz}.mitoALL.filtered.DP.table)
for ((i=3; i<=${ncol}; i +=2)); do
  cut -f $i,$((i+1)) ${VCF%.vcf.gz}.mitoALL.filtered.DP.table |\
  awk '$1 != "./." {print $2}' > $i.mito.DP;
  done

ncol=$(awk 'NR==1{print NF}' ${VCF%.vcf.gz}.nuclearALL.filtered.DP.table)
for ((i=3; i<=${ncol}; i +=2)); do
  cut -f $i,$((i+1)) ${VCF%.vcf.gz}.nuclearALL.filtered.DP.table |\
  awk '$1 != "./." {print $2}' > $i.nuclear.DP;
  done
```

```R
nameList <- c()
for (i in 3:147) { # 21 - odd number for 10 samples
  if (i %% 2==1) nameList <- append(nameList,paste0(i,".mito.DP"))
}

qlist <- matrix(nrow = 73, ncol = 3) # define number of samples (10 samples here)
qlist <- data.frame(qlist, row.names=nameList)
colnames(qlist)<-c('5%', '10%', '99%')

for (i in 1:73) {
  DP <- read.table(nameList[i], header = T)
  qlist[i,] <- quantile(DP[,1], c(.05, .1, .99), na.rm=T)
  d <- density(DP[,1], from=0, to=100, bw=1, na.rm=T)
  ggplot(DP,aes(x=DP[,1])) +
       geom_density() +
       geom_vline(xintercept=c(qlist[i,1],qlist[i,3]), col='red', lwd=1) +
       theme_bw() +
       labs(title = paste0("Sample: ",colnames(DP)), x= "Coverage")
       ggsave(paste0(colnames(DP),"mitoALL.filtered.DP.png"))
}

```
- example of depth plot - these plots have been made for all samples.
![GVCFall.DP.png](../04_analysis/AN_DNK_COG_EN_0012.DPmitoALL.filtered.DP.png)



### Calculating MAD for admix plots
- Guillaume Salle used a MAD approach, calculating admixture for the population, leaving one chromosome out a time.
- Will try this, using scaffolds longer than 1 Mb to make it manageable.

- NOTE: I did not end up using the following code, hence hashed out, but I am keeping it for future reference

```bash
mkdir ~/lustre118_link/trichuris_trichiura/05_ANALYSIS/ADMIXTURE/VALIDATION

cat ../../01_REF/trichuris_trichiura.fa.fai | cut -f1 | grep -v "MITO" > chromosomes.list

# merge the data from individual chromosomes into a single dataset
while read CHROMOSOME; do
     cat $(ls -1 CHROMOSOMES_PL/*BEAGLE.PL | head -n1 | sed 's/CHROMOSOMES_PL\///g') | head -n1 > VALIDATION/no_${CHROMOSOME}.merged.validation.PL;
     cat CHROMOSOMES_PL/*BEAGLE.PL | grep -v "marker" | grep -v "${CHROMOSOME}" >> VALIDATION/no_${CHROMOSOME}.merged.validation.PL;
done < chromosomes.list

# run admixture for multiple values of K

while read CHROMOSOME; do
     for i in 2 3 4 5 6 7 8 9 10; do
     bsub.py --threads 4 2 NGS_admix_multiK "NGSadmix -likes VALIDATION/no_${CHROMOSOME}.merged.validation.PL -K ${i} -P 4 -o VALIDATION/no_${CHROMOSOME}_k_${i}_out -minMaf 0.05" ;
     done;
done  < chromosomes.list

# remove unneeded files
rm *log *gz
```

```R
library(tidyverse)

rep="./"
prefix="no"
k=seq(2,10,1)
chrl=c('Trichuris_trichiura_3_004','Trichuris_trichiura_1_001','Trichuris_trichiura_3_001','Trichuris_trichiura_3_002','Trichuris_trichiura_1_002','Trichuris_trichiura_3_003','Trichuris_trichiura_1_003','Trichuris_trichiura_00_001','Trichuris_trichiura_1_004','Trichuris_trichiura_1_005','Trichuris_trichiura_3_004')
cv1 = array(-1,length(k))
cv2 = array(-1,length(k))
for(K in k){
  for(chr in chrl){
    infile=paste(rep,prefix,"_",chr,"_k_",K,"_out.qopt",sep="")

    # Admixture
    assign(paste0('admix.',chr),t(as.matrix(read.table(infile))))
  }
  # Med absolute deviation
  cv1[K-1] <- mad(c(admix.Trichuris_trichiura_3_004,admix.Trichuris_trichiura_1_001,admix.Trichuris_trichiura_3_001,admix.Trichuris_trichiura_3_002,admix.Trichuris_trichiura_1_002,admix.Trichuris_trichiura_3_003,admix.Trichuris_trichiura_1_003,admix.Trichuris_trichiura_00_001,admix.Trichuris_trichiura_1_004,admix.Trichuris_trichiura_1_005,admix.Trichuris_trichiura_3_004))
  # Jackniffing variance
  cv2[K-1] <- var(c(admix.Trichuris_trichiura_3_004,admix.Trichuris_trichiura_1_001,admix.Trichuris_trichiura_3_001,admix.Trichuris_trichiura_3_002,admix.Trichuris_trichiura_1_002,admix.Trichuris_trichiura_3_003,admix.Trichuris_trichiura_1_003,admix.Trichuris_trichiura_00_001,admix.Trichuris_trichiura_1_004,admix.Trichuris_trichiura_1_005,admix.Trichuris_trichiura_3_004))
  rm(admix.GW_ld.CV1,admix.GW_ld.CV2,admix.GW_ld.CV3,admix.GW_ld.CV4,admix.GW_ld.CV5)
}
```
- ran Guillaumes approach dropping on chromosome at a time, which was ok.
- also reran to generate using different seeds on whole data, and so this is the code to generate the MAD on those data
```R
rep="./"
k=seq(2,10,1)
s=seq(1,5,1)
cv1 = array(-1,length(k))
cv2 = array(-1,length(k))
for(K in k){
  for(S in s){
    infile=paste(rep,"k_",K,"_s_",S,"_out.qopt",sep="")

    # Admixture
    assign(paste0('admix.',S),t(as.matrix(read.table(infile))))
  }
  # Med absolute deviation
  cv1[K-1] <- mad(c(admix.1,admix.2,admix.3,admix.4,admix.5))
  # Jackniffing variance
}

# make a plot
ggplot() +
     geom_point(aes(k,log10(cv1))) +
     geom_line(aes(k,log10(cv1))) +
     theme_bw() +
     labs(x = "Ancestral populations (k)", y = "Median absolute deviation")

# save it
ggsave("plot_admixture_validation_MAD.png")
ggsave("plot_admixture_validation_MAD.pdf", height = 5, width = 5, useDingbats=FALSE)

```
![](../04_analysis/plot_admixture_validation_MAD.png)
```
