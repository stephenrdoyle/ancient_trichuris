# Population genomics of modern and ancient Trichuris trichiura
## Background


1. Reference
2. Processing raw reads
     - trimming
     - kraken check
3. mapping



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


### code
```shell
cd ${WORKING_DIR}/02_RAW

#!/bin/bash
# adaptor remove PE - modern samples
OLD_NAME=${1}
NEW_NAME=${2}

/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/adapterremoval/bin/AdapterRemoval \
--file1 ${OLD_NAME}_R1.merged.fastq.gz \
--file2 ${OLD_NAME}_R2.merged.fastq.gz \
--basename ${NEW_NAME}_PE \
--trimns --trimqualities --collapse --threads 4


#!/bin/bash
# single end - ancient samples
OLD_NAME=${1}
NEW_NAME=${2}

/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/adapterremoval/bin/AdapterRemoval \
--file1 ${OLD_NAME}_R1.merged.fastq.gz \
--basename ${NEW_NAME}_SE \
--trimns --trimqualities --threads 4
```

### run the trimming
```shell
# main samples
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 adapter_remove_modern "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_PE.sh ${OLD_NAME} ${NEW_NAME}"; done < ${WORKING_DIR}/modern.sample_list
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 adapter_remove_ancient "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_SE.sh ${OLD_NAME} ${NEW_NAME}"; done < ${WORKING_DIR}/ancient.sample_list

# other samples
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 adapter_remove_others_PE "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_PE.sh ${OLD_NAME} ${NEW_NAME}"; done < ${WORKING_DIR}/others_PE.sample_list
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 adapter_remove_others_SE "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_SE.sh ${OLD_NAME} ${NEW_NAME}"; done < ${WORKING_DIR}/others_SE.sample_list


```

---


## mapping
Need to map them a little differently. Below are two mapping scripts for each approach.

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


# script for mapping modern samples
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
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 mapping_ancient "${WORKING_DIR}/00_SCRIPTS/run_map_ancient_SE.sh ${OLD_NAME} ${NEW_NAME}" ; done < ${WORKING_DIR}/ancient.sample_list
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 mapping_otherSE "${WORKING_DIR}/00_SCRIPTS/run_map_ancient_SE.sh ${OLD_NAME} ${NEW_NAME}" ; done < ${WORKING_DIR}/others_SE.sample_list

mkdir MAPPED_OTHER
mv OTHER_M* MAPPED_OTHER

mkdir MAPPED_CONTROL
mv CONTROL_* MAPPED_CONTROL


multiqc *flagstat --title mapping
```
[Mapping multiQC report](../04_analysis/mapping_multiqc_report.html)




<!-- ```shell
# generate bamstats
for i in *.bam; do
     bsub.py --queue small --threads 4 2 stats "samtools stats -r ${WORKING_DIR}/01_REF/trichuris_trichiura.fa --threads 4 ${i} \> ${i}.stats" ; done
``` -->


### Kraken
The mapping shows that there is variable mappig rates, and that for some samples there is very poor mapping. This is preticularly the case for the ancient samples, which is to be expected to a degree, given they are both old and collected from the environment. Kraken might give some insight into this, given they might be heavily contaminated with bacteria etc.

```bash
# run kraken on the modern PE trimmed reads
while read OLD_NAME NEW_NAME; do
     bsub.py 10 kraken2 "kraken2 --db /lustre/scratch118/infgen/pathogen/pathpipe/kraken/minikraken_20190423/minikraken2_v1_8GB
     --report ${WORKING_DIR}/02_RAW/${NEW_NAME}.kraken2report
     --paired ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair1.truncated
     ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair2.truncated";
done < ${WORKING_DIR}/modern.sample_list

# run kraken on the ancient SE trimmed reads
while read OLD_NAME NEW_NAME; do
     bsub.py 10 kraken2_SE "kraken2 --db /lustre/scratch118/infgen/pathogen/pathpipe/kraken/minikraken_20190423/minikraken2_v1_8GB
     --report ${WORKING_DIR}/02_RAW/${NEW_NAME}.kraken2report
     ${WORKING_DIR}/02_RAW/${NEW_NAME}_SE.truncated";
done < ${WORKING_DIR}/ancient.sample_list

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




---
