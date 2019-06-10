
## get the raw data
```shell
for i in *gz; do echo ${i%_???.fastq.gz}; done | sort | uniq  |  while read NAME; do cat ${NAME}* > ${NAME}.merged.fastq.gz; done &
```

## trimming
```shell
# PE - modern samples
while read OLD_NAME NEW_NAME; do \
/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/adapterremoval/bin/AdapterRemoval \
--file1 ${OLD_NAME}_R1.merged.fastq.gz \
--file2 ${OLD_NAME}_R2.merged.fastq.gz \
--basename ${NEW_NAME}_PE \
--trimns --trimqualities --collapse --threads 4; done < modern.list



# single end - ancient samples
while read OLD_NAME NEW_NAME; do \
/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/adapterremoval/bin/AdapterRemoval \
--file1 ${OLD_NAME}_R1.merged.fastq.gz \
--basename ${NEW_NAME}_SE \
--trimns --trimqualities --threads 4; done < ancient.list
```
```shell
bsub.py --threads 4 20 adaptor_remove "./adapterremove_2.sh"
```




## script for mapping ancient samples
```shell
#!/bin/bash
# map SE reads
OLD_NAME=${1}
NEW_NAME=${2}

bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/00_REF/trichuris_trichiura.fa ${NEW_NAME}_SE.truncated | samtools view --threads 4 -b - | samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.SE.bam - ;
java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=${NEW_NAME}.tmp.sort.SE.bam OUTPUT=${NEW_NAME}.tmp2.bam METRICS_FILE=${NEW_NAME}.tmp.metrics TMP_DIR=$PWD/${NEW_NAME}.tmp;
samtools flagstat ${NEW_NAME}.tmp2.bam > ${NEW_NAME}.flagstat;
samtools view --threads 4 -F 4 -b -o ${NEW_NAME}.bam ${NEW_NAME}.tmp2.bam;
samtools index -b ${NEW_NAME}.bam;
rm -r ${NEW_NAME}.*tmp*



# script for mapping modern samples
#!/bin/bash
# map PE reads
OLD_NAME=${1}
NEW_NAME=${2}

bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/00_REF/trichuris_trichiura.fa ${NEW_NAME}_PE.pair1.truncated ${NEW_NAME}_PE.pair2.truncated | samtools view --threads 4 -b - | samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.PE.bam - ;
#merge SE reads
cat ${NEW_NAME}_PE.singleton.truncated ${NEW_NAME}_PE.collapsed.truncated ${NEW_NAME}_PE.collapsed > ${NEW_NAME}_SE.tmp.fastq ;
# map SE reads
bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/00_REF/trichuris_trichiura.fa ${NEW_NAME}_SE.tmp.fastq | samtools view --threads 4 -b - | samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.SE.bam - ;
samtools merge ${NEW_NAME}.tmp.bam ${NEW_NAME}.tmp.sort.PE.bam ${NEW_NAME}.tmp.sort.SE.bam;
java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=${NEW_NAME}.tmp.bam OUTPUT=${NEW_NAME}.tmp2.bam METRICS_FILE=${NEW_NAME}.tmp.metrics TMP_DIR=$PWD/${NEW_NAME}.tmp;
samtools flagstat ${NEW_NAME}.tmp2.bam > ${NEW_NAME}.flagstat;
samtools view --threads 4 -F 4 -b -o ${NEW_NAME}.bam ${NEW_NAME}.tmp2.bam;
samtools index -b ${NEW_NAME}.bam;
rm -r ${NEW_NAME}.*tmp*
```


```shell
# run the mapping jobs
while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 modern "./run_map_modern.sh ${OLD_NAME} ${NEW_NAME}" ; done < modern.list

while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 ancient "./run_map_ancient.sh ${OLD_NAME} ${NEW_NAME}" ; done < ancient.list
```
