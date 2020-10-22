# Population genomics of modern and ancient Trichuris trichiura

## Contents

1. Background
2. Reference
3. Processing raw reads
     - trimming
4. mapping
     - kraken check of raw reads post mapping
5. QC
     - deamination damage
     - trim bams to correct for damage



***
## Background


### Sampling (from martins slides - need to check this against samples i have)
- Modern genomes  from human hosts (49)
     - Human worm isolates (35)
          - Uganda (12), China (7), Ecuador (8), Honduras (8)
     - Human egg isolates (14)
          - Cameroon (5), Ethiopia (4), Tanzania (5)

- Modern genomes  from non-human hosts (7)
     - Worm isolates (7)
          - Baboon – Denmark (2), Colobus - Spain (2), Leaf-Monkey – China (3)

- Ancient samples (19 samples from 10 sites)
     - 4-500 BC: Qala’at al-Bahrain (1) –  in analysis - **can't find this data**
     - 300 BC: Tollundmanden, DK (1) – in analysis - **can't find this data**
     - 1000 AD: Viborg, DK (2) – 2x Mito + 0,01x/0,16x NG **can't find this data**
     - 1350 AD : Kampen, Netherlands (2) – 2x Mito + 1,8x/4,7x NG - yes
     - 1350-1400 AD: Odense, DK (4) – in analysis - yes
     - 1500-1600 AD: Zwolle, Netherlands (1) – Mito + 2,3x NG - yes
     - 1550-1580 AD: Vilnius, Lithuania (1) – in analysis - yes
     - 1600-1800 AD: Gl. Strand, CPH, DK (4) – in analysis - yes
     - 1680 AD: Kultorvet, CPH, DK (2) – in analysis - yes
     - 1700 AD: Adelgade, CPH, DK (1) – in analysis - yes


### Overarching questions
- describe new genome?
- how similar are worms from humans and animals?
- how similar are modern and ancient samples?
- global dispersal timing?
- evidence of selection around beta-tubulin?>
     - might apply to some modern samples


### Downstream analyses
- PCA
     - whole genomes modern
     - mtDNA genomes all
- faststructure
- treemix
- MSMC
- genome-wide patterns of genetic diversity
     - Pi
     - Tajimas D
     - Fst





***
## working directory
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura
WORKING_DIR=/nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura

```


## Project setup
```shell
mkdir 00_SCRIPTS 01_REF 02_RAW 03_MAPPING
```
---

## Reference
- the reference is the unpublished Trichuris trichiura assembly
- genome stats
     - sum = 80573711, n = 113, ave = 713041.69, largest = 29164577
     - N50 = 11299416, n = 2
     - N60 = 9167782, n = 3
     - N70 = 5100676, n = 5
     - N80 = 2017282, n = 7
     - N90 = 643749, n = 14
     - N100 = 1808, n = 113
     - N_count = 250000
     - Gaps = 25

```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/01_REF

# make a generic bwa index for mapping later on
bwa index trichuris_trichiura.fa

# make a generic dict file for SNP calling later on
samtools dict trichuris_trichiura.fa > trichuris_trichiura.dict
```
---

## get the raw data
- I received the raw sequnecing data from Peter on a hard drive, and copied all to the lustre environment as is into the directory 02_RAW
- want to collect all of the fastq files into one place, and then begin processing them

```shell
cd ${WORKING_DIR}/02_RAW

for i in *gz; do echo ${i%_???.fastq.gz}; done | sort | uniq  |  while read NAME; do cat ${NAME}* > ${NAME}.merged.fastq.gz; done &

# sample used in the Foth 2014 genome paper
pf data --type lane --id 6929_8 -l ./ --filetype fastq

```
- this generated a total of 130 R1 files and 60 R2 files
- the difference in number reflects the fact that the "modern" samples have been sequenced using PE reads whereas the "ancient" samples have been sequenced using SE reads. Will need to treat these a little differently later on.
- there is also a difference in the number of samples with useable data (based on some sample informaiton spreadsheets from MS / PJ), and the total number
     - some samples are controls
     - some samples presumably didnt work, and so we excluded from the preliminary analyses
     - some samples are duplicates, and so will need to be merged at some point too.
---

## metadata
- I have begun to curate the metadata, to get a better idea of what data is actually available, and to ensure that I have consistent informaiton for all samples.
- The existing naming scheme is not very informative to me at least, and so will rename everything with a simple format
- GoogleSheet:  https://docs.google.com/spreadsheets/d/1PiapiaZZw0g0i3lN0feXxqEVupvaAuOm0IViXnZPFqk/edit?usp=sharing

---
## trimming
- trimming using AdapterRemoval, which seems to be used for a few different ancient DNA projects. I think it is because it trims and merges, which generally improves the mapping scores off some poor qual end of reads.
- Tool: https://buildmedia.readthedocs.org/media/pdf/adapterremoval/latest/adapterremoval.pdf

```bash
cd ${WORKING_DIR}/02_RAW

### run the trimming
```shell
# main samples
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 adapter_remove_modern "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_PE.sh ${OLD_NAME} ${NEW_NAME}"; done < ${WORKING_DIR}/modern.sample_list
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 adapter_remove_ancient "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_SE.sh ${OLD_NAME} ${NEW_NAME}"; done < ${WORKING_DIR}/ancient.sample_list

# other samples
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 adapter_remove_others_PE "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_PE.sh ${OLD_NAME} ${NEW_NAME}"; done < ${WORKING_DIR}/others_PE.sample_list
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 adapter_remove_others_SE "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_SE.sh ${OLD_NAME} ${NEW_NAME}"; done < ${WORKING_DIR}/others_SE.sample_list

```
where "run_adapter_remove_PE.sh" is:
```bash
#!/bin/bash
# adaptor remove PE - modern samples
OLD_NAME=${1}
NEW_NAME=${2}

/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/adapterremoval/bin/AdapterRemoval \
--file1 ${OLD_NAME}_R1.merged.fastq.gz \
--file2 ${OLD_NAME}_R2.merged.fastq.gz \
--basename ${NEW_NAME}_PE \
--trimns --trimqualities --collapse --threads 4

```

and where "run_adapter_remove_SE.sh" is:
```bash
#!/bin/bash
# single end - ancient samples
OLD_NAME=${1}
NEW_NAME=${2}

/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/adapterremoval/bin/AdapterRemoval \
--file1 ${OLD_NAME}_R1.merged.fastq.gz \
--basename ${NEW_NAME}_SE \
--trimns --trimqualities --threads 4
```



### FIX: merge duplicate ancient read sets
These seem to have been a single sample, extracted twice (or in two ways, perhaps two washes of a column) and sequenced individually. Given the low coverage, these should be merged (I think).
```bash
# merge the duplicated set, renaming with a unique name that refers to both the originals, and then remove the originals
cat AN_DNK_COG_EN_001_SE.truncated AN_DNK_COG_EN_002_SE.truncated > AN_DNK_COG_EN_0012_SE.truncated; rm AN_DNK_COG_EN_001_SE.truncated AN_DNK_COG_EN_002_SE.truncated
cat AN_DNK_COG_EN_003_SE.truncated AN_DNK_COG_EN_004_SE.truncated > AN_DNK_COG_EN_0034_SE.truncated; rm AN_DNK_COG_EN_003_SE.truncated AN_DNK_COG_EN_004_SE.truncated
cat AN_DNK_COG_EN_005_SE.truncated AN_DNK_COG_EN_006_SE.truncated > AN_DNK_COG_EN_0056_SE.truncated; rm AN_DNK_COG_EN_005_SE.truncated AN_DNK_COG_EN_006_SE.truncated
cat AN_DNK_COG_EN_007_SE.truncated AN_DNK_COG_EN_008_SE.truncated > AN_DNK_COG_EN_0078_SE.truncated; rm AN_DNK_COG_EN_007_SE.truncated AN_DNK_COG_EN_008_SE.truncated
cat AN_DNK_COK_EN_001_SE.truncated AN_DNK_COK_EN_002_SE.truncated > AN_DNK_COK_EN_0012_SE.truncated; rm AN_DNK_COK_EN_001_SE.truncated AN_DNK_COK_EN_002_SE.truncated
cat AN_DNK_COK_EN_003_SE.truncated AN_DNK_COK_EN_004_SE.truncated > AN_DNK_COK_EN_0034_SE.truncated; rm AN_DNK_COK_EN_003_SE.truncated AN_DNK_COK_EN_004_SE.truncated
cat AN_DNK_OBM_EN_001_SE.truncated AN_DNK_OBM_EN_002_SE.truncated > AN_DNK_OBM_EN_0012_SE.truncated; rm AN_DNK_OBM_EN_001_SE.truncated AN_DNK_OBM_EN_002_SE.truncated
cat AN_DNK_OBM_EN_003_SE.truncated AN_DNK_OBM_EN_004_SE.truncated > AN_DNK_OBM_EN_0034_SE.truncated; rm AN_DNK_OBM_EN_003_SE.truncated AN_DNK_OBM_EN_004_SE.truncated
cat AN_DNK_OBM_EN_005_SE.truncated AN_DNK_OBM_EN_006_SE.truncated > AN_DNK_OBM_EN_0056_SE.truncated; rm AN_DNK_OBM_EN_005_SE.truncated AN_DNK_OBM_EN_006_SE.truncated
cat AN_DNK_OBM_EN_007_SE.truncated AN_DNK_OBM_EN_008_SE.truncated > AN_DNK_OBM_EN_0078_SE.truncated; rm AN_DNK_OBM_EN_007_SE.truncated AN_DNK_OBM_EN_008_SE.truncated
cat AN_DNK_OBM_EN_009_SE.truncated AN_DNK_OBM_EN_010_SE.truncated > AN_DNK_OBM_EN_0910_SE.truncated; rm AN_DNK_OBM_EN_009_SE.truncated AN_DNK_OBM_EN_010_SE.truncated
cat AN_LTU_VIL_EN_001_SE.truncated AN_LTU_VIL_EN_002_SE.truncated > AN_LTU_VIL_EN_0012_SE.truncated; rm AN_LTU_VIL_EN_001_SE.truncated AN_LTU_VIL_EN_002_SE.truncated
cat AN_NLD_KAM_EN_001_SE.truncated AN_NLD_KAM_EN_002_SE.truncated > AN_NLD_KAM_EN_0012_SE.truncated; rm AN_NLD_KAM_EN_001_SE.truncated AN_NLD_KAM_EN_002_SE.truncated
cat AN_NLD_KAM_EN_003_SE.truncated AN_NLD_KAM_EN_004_SE.truncated > AN_NLD_KAM_EN_0034_SE.truncated; rm AN_NLD_KAM_EN_003_SE.truncated AN_NLD_KAM_EN_004_SE.truncated
cat AN_NLD_ZWO_EN_001_SE.truncated AN_NLD_ZWO_NA_002_SE.truncated > AN_NLD_ZWO_EN_0012_SE.truncated; rm AN_NLD_ZWO_EN_001_SE.truncated AN_NLD_ZWO_NA_002_SE.truncated

# make a new sample list to work from using new names
ls -1 AN* | cut -c-18 > ../ancient.sample_list_v2
```


```bash
# clean up
rm *discarded *settings

```




---


## mapping
Need to map them a little differently. Below are two mapping scripts for each approach.

### script for mapping the ancient samples - these are all single end (SE) reads
```shell
#!/bin/bash
# map SE reads
OLD_NAME=${1}
NEW_NAME=${2}

# map the SE reads
bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M ${WORKING_DIR}/01_REF/trichuris_trichiura.fa ${WORKING_DIR}/02_RAW/${NEW_NAME}_SE.truncated |\
     samtools view --threads 4 -b - |\
     samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.SE.bam - ;

# mark duplicate reads
java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=${NEW_NAME}.tmp.sort.SE.bam OUTPUT=${NEW_NAME}.tmp2.bam METRICS_FILE=${NEW_NAME}.tmp.metrics TMP_DIR=$PWD/${NEW_NAME}.tmp;

# generate mapping stats from unfilted bam
samtools flagstat ${NEW_NAME}.tmp2.bam > ${NEW_NAME}.flagstat;

# filter bam to only contain mapped reads
samtools view --threads 4 -F 4 -b -o ${NEW_NAME}.bam ${NEW_NAME}.tmp2.bam;

# index the filtered bam
samtools index -b ${NEW_NAME}.bam;

# clean up unnecessary tmp files
rm -r ${NEW_NAME}.*tmp*
```


### script for mapping modern samples - these are all paired-end (PE) reads
```shell
#!/bin/bash
# map PE reads
OLD_NAME=${1}
NEW_NAME=${2}

# map the PE reads
bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M ${WORKING_DIR}/01_REF/trichuris_trichiura.fa ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair1.truncated ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair2.truncated |\
     samtools view --threads 4 -b - |\
     samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.PE.bam - ;

#merge SE reads that are either merged PE reads, or singletons after adapter filtering
cat ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.singleton.truncated ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.collapsed.truncated ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.collapsed > ${NEW_NAME}_SE.tmp.fastq ;

# map SE reads
bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M ${WORKING_DIR}/01_REF/trichuris_trichiura.fa ${NEW_NAME}_SE.tmp.fastq | \
     samtools view --threads 4 -b - |\
     samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.SE.bam - ;

# merge PE and SE bams
samtools merge ${NEW_NAME}.tmp.bam ${NEW_NAME}.tmp.sort.PE.bam ${NEW_NAME}.tmp.sort.SE.bam;

# mark duplicate reads
java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=${NEW_NAME}.tmp.bam OUTPUT=${NEW_NAME}.tmp2.bam METRICS_FILE=${NEW_NAME}.tmp.metrics TMP_DIR=$PWD/${NEW_NAME}.tmp;

# generate mapping stats from unfilted bam
samtools flagstat ${NEW_NAME}.tmp2.bam > ${NEW_NAME}.flagstat;

# filter bam to only contain mapped reads
samtools view --threads 4 -F 4 -b -o ${NEW_NAME}.bam ${NEW_NAME}.tmp2.bam;

# index the filtered bam
samtools index -b ${NEW_NAME}.bam;

# clean up unnecessary tmp files
rm -r ${NEW_NAME}.*tmp*
```


```shell
cd ${WORKING_DIR}/03_MAPPING

# run the mapping jobs for modern and ancient samples
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 mapping_modern "${WORKING_DIR}/00_SCRIPTS/run_map_modern_PE.sh ${OLD_NAME} ${NEW_NAME}" ; done < ${WORKING_DIR}/modern.sample_list
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 mapping_otherPE "${WORKING_DIR}/00_SCRIPTS/run_map_modern_PE.sh ${OLD_NAME} ${NEW_NAME}" ; done < ${WORKING_DIR}/others_PE.sample_list

# run the mapping jobs for the control and other samples
#while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 mapping_ancient "${WORKING_DIR}/00_SCRIPTS/run_map_ancient_SE.sh ${OLD_NAME} ${NEW_NAME}" ; done < ${WORKING_DIR}/ancient.sample_list_v2
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 mapping_otherSE "${WORKING_DIR}/00_SCRIPTS/run_map_ancient_SE.sh ${OLD_NAME} ${NEW_NAME}" ; done < ${WORKING_DIR}/others_SE.sample_list

mkdir MAPPED_OTHER
mv OTHER_M* MAPPED_OTHER

mkdir MAPPED_CONTROL
mv CONTROL_* MAPPED_CONTROL


# rerun of ancient samples due to name change
while read NEW_NAME; do bsub.py --threads 4 20 mapping_ancient "${WORKING_DIR}/00_SCRIPTS/run_map_ancient_SE.sh NULL ${NEW_NAME}" ; done < ${WORKING_DIR}/ancient.sample_list_v2



multiqc *flagstat --title mapping
```
[Mapping MultiQC report](../04_analysis/mapping_multiqc_report.html)




<!-- ```shell
# generate bamstats
for i in *.bam; do
     bsub.py --queue small --threads 4 2 stats "samtools stats -r ${WORKING_DIR}/01_REF/trichuris_trichiura.fa --threads 4 ${i} \> ${i}.stats" ; done
``` -->


### Kraken of trimmed reads post mapping
The mapping shows that there is variable mapping rates, and that for some samples there is very poor mapping. This is particularly the case for the ancient samples, which is to be expected to a degree, given they are both old and collected from the environment. Kraken might give some insight into this, given they might be heavily contaminated with bacteria etc.

```bash
module load kraken2/2.0.8_beta

# run kraken on the modern PE trimmed reads
while read OLD_NAME NEW_NAME; do
     bsub.py 10 kraken2 "kraken2 --db /lustre/scratch118/infgen/pathogen/pathpipe/kraken/minikraken_20190423/minikraken2_v1_8GB \
     --report ${WORKING_DIR}/02_RAW/${NEW_NAME}.kraken2report \
     --paired ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair1.truncated \
     ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair2.truncated";
done < ${WORKING_DIR}/modern.sample_list

# run kraken on the ancient SE trimmed reads
# while read OLD_NAME NEW_NAME; do
#      bsub.py 10 kraken2_SE "kraken2 --db /lustre/scratch118/infgen/pathogen/pathpipe/kraken/minikraken_20190423/minikraken2_v1_8GB
#      --report ${WORKING_DIR}/02_RAW/${NEW_NAME}.kraken2report
#      ${WORKING_DIR}/02_RAW/${NEW_NAME}_SE.truncated";
# done < ${WORKING_DIR}/ancient.sample_list

while read NEW_NAME; do
     bsub.py 5 kraken2_SE "kraken2 --db /lustre/scratch118/infgen/pathogen/pathpipe/kraken/minikraken_20190423/minikraken2_v1_8GB \
     --report ${WORKING_DIR}/02_RAW/${NEW_NAME}.kraken2report \
     ${WORKING_DIR}/02_RAW/${NEW_NAME}_SE.truncated";
done < ${WORKING_DIR}/ancient.sample_list_v2

# once the kraken runs have completed, run multiqc .
multiqc *kraken2report --title kraken
```
[Kraken MultiQC report](../04_analysis/kraken_multiqc_report.html)
- the output shows most samples have a small degree of contamination based on hits in the kraken database
- non have a lot of contamination, whcih is slightly surprising
- this alone doesnt explain the mismapping, althoguh, it simply may mean that the putative contaminant is not present in the kraken databased
- could try
     - blasting some sequences
     - mapping to other reference - are there other nematodes in these environmental ancient samples?






## Damage
```bash
# To view deamination-derived damage patterns in a simple table, without separating CpG sites
#samtools view AN_DNK_COG_EN_002.bam | python pmdtools.0.60.py --deamination

# modern samples
while read -r OLD_NAME NEW_NAME; do
# To compute deamination-derived damage patterns separating CpG and non-CpG sites
samtools view ${NEW_NAME}.bam | head -n 10000 | pmdtools --platypus --requirebaseq 30 > PMD_temp.txt ;

R CMD BATCH ${WORKING_DIR}/00_SCRIPTS/plotPMD.R ;

mv deamination_plot.png ${NEW_NAME}.deamination_plot.png ; done < ${WORKING_DIR}/modern.sample_list

# ancient samples
while read -r NEW_NAME; do
# To compute deamination-derived damage patterns separating CpG and non-CpG sites
samtools view ${NEW_NAME}.bam | head -n 10000 | pmdtools --platypus --requirebaseq 30 > PMD_temp.txt ;

R CMD BATCH ${WORKING_DIR}/00_SCRIPTS/plotPMD.R ;

mv deamination_plot.png ${NEW_NAME}.deamination_plot.png ; done < ancient.sample_list_v2



# once finished, move all the plots into a new directory
mkdir ${WORKING_DIR}/04_ANALYSES/deamination && mv ${WORKING_DIR}/03_MAPPING/*deamination_plot.png ${WORKING_DIR}/04_ANALYSES/deamination
```

where "plotPMD.R" is:
```R
# script to make deamination frequency plots from PMD data

# load libraries
require(reshape2)
require(patchwork)
require(ggplot2)
require(dplyr)
require(tidyverse)

# read data from PMD, example "samtools view AN_DNK_COG_EN_002.bam | head -n 100000 | pmdtools --platypus --requirebaseq 30 > PMD_temp.txt"
data <- read.table("PMD_temp.txt", header=T)

data2 <- melt(data,id.vars="z")

# split into 5' and 3' datasets, and add some variables for plotting
data_5 <- data2 %>%
     filter(str_detect(variable, "5")) %>%
     mutate(base_substitution = if_else(str_detect(variable,"GA"), "G-to-A", if_else(str_detect(variable,"CT"), "C-to-T", "other"))) %>%
     mutate(CpG_state = if_else(str_detect(variable,"CpG"), "CpG", "non-CpG"))

data_3 <- data2 %>%
     filter(str_detect(variable, "3")) %>%
     mutate(base_substitution = if_else(str_detect(variable,"GA"), "G-to-A", if_else(str_detect(variable,"CT"), "C-to-T", "other"))) %>%
     mutate(CpG_state = if_else(str_detect(variable,"CpG"), "CpG", "non-CpG"))

# make plots
plot_5 <- ggplot(data_5, aes(z, value, colour = base_substitution, group = variable, linetype = CpG_state, size = base_substitution)) +
     geom_line() +
     ylim(0, 0.2) +
     theme_bw() +
     scale_colour_manual(values = c("other" = "black", "C-to-T" = "red" , "G-to-A" = "blue"))+
     scale_size_manual(values=c("other" = 0.5, "C-to-T" = 1 , "G-to-A" = 1))+
     scale_linetype_manual(values = c("CpG" = "dashed", "non-CpG" = "solid")) +
     labs(x="Distance from 5' end of sequence read", y="Mismatch frequency")

plot_3 <- ggplot(data_3, aes(z, value, colour = base_substitution, group = variable, linetype = CpG_state, size = base_substitution)) +
     geom_line() +
     ylim(0, 0.2) +
     theme_bw() +
     scale_colour_manual(values = c("other" = "black", "C-to-T" = "red" , "G-to-A" = "blue"))+
     scale_size_manual(values=c("other" = 0.5, "C-to-T" = 1 , "G-to-A" = 1))+
     scale_linetype_manual(values = c("CpG" = "dashed", "non-CpG" = "solid")) +
     labs(x="Distance from 3' end of sequence read", y="Mismatch frequency")

# bring it together
plot_5 + plot_3 + plot_layout(guides = "collect")

# save it
ggsave("deamination_plot.png", height=5, width=10)

```

- Example of a modern sample (MN_CHN_GUA_HS_001)
![MN_CHN_GUA_HS_001.deamination_plot](../04_analysis/MN_CHN_GUA_HS_001.deamination_plot.png)
- Example of a ancient sample (AN_DNK_COG_EN_002)  
![AN_DNK_COG_EN_002.deamination_plot](../04_analysis/AN_DNK_COG_EN_002.deamination_plot.png)

- clearly a CT bias in the first two bases of the ancient sample, that doesnt seem to be present in the modern sample
     - not present in all samples I dont think, but in quite a few that I have checked.
- simplest solution is to remove the first two bases from all reads before moving forward




### trim bases from reads in bam
```bash
while read -r OLD_NAME NEW_NAME; do
     bsub.py 10 --threads 4 trim_bams "${WORKING_DIR}/00_SCRIPTS/run_trimreads_in_bam.sh ${NEW_NAME}";
done < ${WORKING_DIR}/modern.sample_list

while read -r NEW_NAME; do
     bsub.py 10 --threads 4 trim_bams "${WORKING_DIR}/00_SCRIPTS/run_trimreads_in_bam.sh ${NEW_NAME}";
done < ${WORKING_DIR}/ancient.sample_list_v2

```

where "run_trimreads_in_bam.sh" is:
```
#!/bin/bash

# trim left and right bases from reads in a bam.

NAME=${1}

bamutils_clip_left=2
bamutils_clip_right=2


bamUtils trimBam ${NAME}.bam tmp.bam -L ${bamutils_clip_left} -R ${bamutils_clip_right}
samtools sort -@ 4 tmp.bam -o ${NAME}.trimmed.bam
samtools index ${NAME}.trimmed.bam
```










## World sampling map
Given it is a "global diversity" study, worth having a world map with sampling sites, distinction between ancient and modern samples, and the fact that some some from humans, animals, and the environment (ancient).
```R
setwd("/nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/05_ANALYSIS/MAP")

# load libraries
library(ggplot2)
library(dplyr)
require(maps)
require(viridis)
library(ggrepel)

# load world data
world_map <- map_data("world")

# load metadata
data <- read.table("map_metadata.txt", header=T, sep="\t")

# make a map
ggplot() +
  geom_polygon(data = world_map, aes(x = world_map$long, y = world_map$lat, group = world_map$group), fill="grey90") +
  geom_point(data = data, aes(x = LONG, y = LAT, colour = SAMPLE_AGE, shape = SAMPLE_LOCATION), size=3) +
  geom_text_repel(data = data, aes(x = LONG, y = LAT, label = paste0(COUNTRY," (",REGION,"); n = ", SAMPLE_N)), size=3) +        
  theme_void() +
  ylim(-55,85) +
  labs(colour="", shape="")

# save it
ggsave("worldmap_samplingsites.png", height=5, width=12)
#ggsave("worldmap_samplingsites.pdf", height=5, width=12)
```
![worldmap_samplingsites](../04_analysis/worldmap_samplingsites.png)

---
