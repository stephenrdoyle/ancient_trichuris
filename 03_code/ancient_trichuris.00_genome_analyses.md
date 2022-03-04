# Genome analyses

## BUSCO


```bash

export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/pkgs/augustus-3.1-0/bin:$PATH"
export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/pkgs/augustus-3.1-0/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/pkgs/augustus-3.1-0/config"

bsub.py --threads 20 --queue long 20 busco_tt_new "busco -i trichuris_trichiura.fa -l /nfs/users/nfs_s/sd21/lustre118_link/databases/busco/metazoa_odb9 -o new_assembly_metazoa --mode genome --long -sp caenorhabditis -f --cpu 20"

bsub.py --threads 20 --queue long 20 busco_tt_old "busco -i trichuris_trichiura.PRJEB535.WBPS12.genomic.fa -l /nfs/users/nfs_s/sd21/lustre118_link/databases/busco/metazoa_odb9 -o old_assembly_metazoa --mode genome --long -sp caenorhabditis -f --cpu 20"

```

- Summarized benchmarking in BUSCO notation for file trichuris_trichiura.PRJEB535.WBPS12.genomic.fa
- BUSCO was run in mode: genome
     - C:81.9%[S:79.8%,D:2.1%],F:1.8%,M:16.3%,n:978
          - 801	Complete BUSCOs (C)
          - 780	Complete and single-copy BUSCOs (S)
          - 21	Complete and duplicated BUSCOs (D)
          - 18	Fragmented BUSCOs (F)
          - 159	Missing BUSCOs (M)
          - 978	Total BUSCO groups searched


- Summarized benchmarking in BUSCO notation for file trichuris_trichiura.fa
- BUSCO was run in mode: genome
     - C:81.4%[S:79.3%,D:2.1%],F:1.7%,M:16.9%,n:978
          - 797	Complete BUSCOs (C)
          - 776	Complete and single-copy BUSCOs (S)
          - 21	Complete and duplicated BUSCOs (D)
          - 17	Fragmented BUSCOs (F)
          - 164	Missing BUSCOs (M)
          - 978	Total BUSCO groups searched







## Liftoff
- need to identify genes in new genome assembly, which does not yet have a stable annotations
- propose to simply lift over gene models from original annotation, simply for the purpose of asking, what genes are present
- doesnt have to be perfect
- using / trying the new "liftoff" tool - https://github.com/agshumate/Liftoff; https://doi.org/10.1093/bioinformatics/btaa1016


```bash

conda activate py37

pip3 install liftoff --user

conda install -c bioconda minimap2


wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/trichuris_trichiura/PRJEB535/trichuris_trichiura.PRJEB535.WBPS15.annotations.gff3.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/trichuris_trichiura/PRJEB535/trichuris_trichiura.PRJEB535.WBPS15.genomic.fa.gz

liftoff -g trichuris_trichiura.PRJEB535.WBPS15.annotations.gff3.gz trichuris_trichiura.fa trichuris_trichiura.PRJEB535.WBPS15.genomic.fa -o newversion.gff -u unmapped_features.txt


wc -l unmapped_features.txt
#> 1199 unmapped_features.txt

awk $3="gene" '{print}' newversion.gff
#> 8451 genes



# to be cross-compatible with Apollo, which has the same genome version but with updated scaffold names (changed by Faye after I had started the analysis), had to run the following on the genome and annotation

perl rename_scaffolds.pl trichuris_trichiura.fa > trichuris_trichiura.renamed.fa

perl rename_scaffolds.pl newversion.gff > newversion.renamed.gff

# where "rename_scaffolds.pl" is:

#!/usr/bin/env perl
use strict;
use warnings;
while (<>){
    s/Trichuris_trichiura/TTRE/;
    s/_00_/_unplaced_/;
    s/_1_/_chr1_/;
    s/_2_001/_chr2/;
    s/_3_/_chr3_/;
    s/_[0]+/_scaffold/;
    print $_;
}

# found that, becasue the original version has been annotated interpro descriptions etc, the liftover has misformatted some of the lines (in the same way that apollo dumps of the Haemonchus annotations were misformatted). To fix, ran the following:

cat newversion.renamed.gff | while read LINE; do
     if [[ $LINE =~ (^TTRE*|\#) ]]; then
          echo -ne "\n${LINE} ";
          else
          echo -n "${LINE} ";
     fi;
done > newversion.renamed.v2.gff3


cat newversion.gff | while read LINE; do
     if [[ $LINE =~ (^Trichuris*|\#) ]]; then
          echo -ne "\n${LINE} ";
          else
          echo -n "${LINE} ";
     fi;
done > newversion.v2.gff3

echo -e "##gff-version 3" > UPDATED_annotation.gff3
cat newversion.renamed.v2.gff3 | sort -k1,1 -k4,4n | sed '/^$/d' >> UPDATED_annotation.gff3

```

- ran really quickly, maybe a few minutes at most
- shows 1199 features in the "unmapped_features.txt" - have not idea how "bad" this is, but possible ok(?) coming from a draft assembly to something that is more complete. Very likely some at least are errors, but possibly not a big deal for my use.
- 8451 genes

- will attempt to map protein sequences back to genome with exonerate

```bash
# extract protein sequences
gffread <(zcat trichuris_trichiura.PRJEB535.WBPS15.annotations.gff3.gz) -g trichuris_trichiura.PRJEB535.WBPS15.genomic.fa -y old_proteins.fa
gffread <(zcat trichuris_trichiura.PRJEB535.WBPS15.annotations.gff3.gz) -g trichuris_trichiura.PRJEB535.WBPS15.genomic.fa -x old_cds.fa

# fix line lengths
fastaq to_fasta -l0 old_proteins.fa tmp; mv tmp old_proteins.fa
fastaq to_fasta -l0 old_cds.fa tmp; mv tmp old_cds.fa

# extract protein sequences for the unmapped genes
cat unmapped_features.txt | sed 's/'gene:'//g' | while read -r GENE; do
     grep -A1 ${GENE} old_proteins.fa >> unmapped_proteins.fa;
done

cat unmapped_features.txt | sed 's/'gene:'//g' | while read -r GENE; do
     grep -A1 ${GENE} old_cds.fa >> unmapped_cds.fa;
done

~sd21/bash_scripts/run_exonerate_splitter trichuris_trichiura.fa unmapped_proteins.fa

cat split_exonerate*out | Exonerate_to_evm_gff3.pl - > merged_exonerate.output

rm split_exonerate*out
rm x*
rm run_split*

```



# run diamond to query species
```bash
# load diamond
module load diamond/0.9.24--ha888412_1

diamond blastp \
--db /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/uniprot_ref_proteomes.dmnd \
--query unmapped_proteins.fa \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
--taxonnodes /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/nodes.dmp \
--tmpdir /dev/shm \
--taxonmap /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/prot.accession2taxid \
--max-target-seqs 1 > unmapped_proteins.diamond.out

diamond blastx \
--db /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/uniprot_ref_proteomes.dmnd \
--query unmapped_cds.fa \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
--taxonnodes /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/nodes.dmp \
--tmpdir /dev/shm \
--taxonmap /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/prot.accession2taxid \
--max-target-seqs 1 > unmapped_cds.diamond.out


cut -f4 unmapped_proteins.diamond.out | sort | uniq -c | sort -k1nr

```

796 of 1000 hits classified as contaminants
204 not in uniprot refseq database

588 83333 Escherichia coli K-12
204  No hit
161 224308     Bacillus subtilis
13 6239   Caenorhabditis elegans
 3 10090  Mus musculus
 3 224911 Bradyrhizobium diazoefficiens
 3 243274 Thermotoga maritima
 2 100226 Streptomyces coelicolor
 2 1111708     Synechocystis sp
 2 122586 Neisseria meningitidis
 2 190304 Fusobacterium nucleatum
 2 208964 Pseudomonas aeruginosa
 2 224324 Aquifex aeolicus
 2 7955   Danio rerio
 1 243090 Rhodopirellula baltica
 1 243230 Deinococcus radiodurans
 1 243232 Methanocaldococcus jannaschii
 1 273057 Saccharolobus solfataricus
 1 44689 Dictyostelium discoideum
 1 515635 Dictyoglomus turgidum
 1 7227 Drosophila melanogaster
 1 83332 Mycobacterium tuberculosis
 1 8364 Xenopus tropicalis
 1 85962 Helicobacter pylori
 1 9913 Bos taurus

 cut -f13 unmapped_cds.diamond.out | sort | uniq -c | sort -k1nr
     598 83333 Escherichia coli K-12
     206
     161 224308  Bacillus subtilis
      13 6239
       3 10090
       3 224911
       3 243274
       2 100226
       2 1111708
       2 122586
       2 190304
       2 208964
       2 224324
       2 7955
       1 243090
       1 243230
       1 243232
       1 273057
       1 44689
       1 515635
       1 7227
       1 83332
       1 8364
       1 85962
       1 9913




## Running interproscan on liftoff annotation and making a list of Go terms per gene

```bash
# make a proteins file from annotation
gffread -y PROTEINS.fa -g trichuris_trichiura.fa liftover_annotation.gff3

# remove stop codons from protein sequences - IPS doesnt like it
sed -e 's/\.//g' PROTEINS.fa > tmp; mv tmp PROTEINS.fa

# load module
module load interproscan/5.39-77.0-W01

# run IPS
farm_interproscan -a PROTEINS.fa -o IPS.output.gff

# lift over GO terms from interproscan to GFF
extract_interproscan_go_terms -i IPS.output.gff -e liftover_annotation.gff3



# extract GO terms from mRNAs in GFF
awk '$3=="mRNA" {print $0}' OFS="\t" liftover_annotation.gff3.go.gff | grep "Ontology_term" | cut -f 9 | cut -f3,4 -d ";" | sed -e 's/;/\t/g' -e 's/Name=//g' -e 's/"//g' -e 's/Ontology_term=//g' -e 's/-.*\t/\t/g' -e 's/,/\t/g' -e 's/\..*\t/\t/g' > annotation_GO_per_gene.txt

# work through columns to split multiple go terms per gene into one term per gene, repeating the gene name if multiple GO terms present
awk '{for(i=2; i<=NF; i++) {print $1,$i}}' OFS="\t" annotation_GO_per_gene.txt > annotation_GO_per_gene_split.txt

```









## mitochondrial Genome annotation
- first, want to map get protein coding gene annotations onto reference assembly, which doesnt have any annotations
- using liftoff () to transfer from
- didnt realise this, but the mito genome in the assembly is quite different to Foth et al.
     - only about 81% similarity
- blastn at NCBI identified KT449826.1 as a much better match
     - this had an annotation, so using that as the liftover
     - downloaded both the fasta and GFF3 from NCBI

```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/01_REF

conda activate liftoff

liftoff -g KT449826.1.gff3 Trichuris_trichiura_MITO.fasta KT449826.1.fa -o Trichuris_trichiura_MITO.gff3

```

# make a mask of non-protein-coding regions
```bash
cat Trichuris_trichiura_MITO.fasta.fai | cut -f1,2 > mito.genome

bedtools complement -i Trichuris_trichiura_MITO.gff3 -g mito.genome > mito.mask.bed

```


```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/05_ANALYSIS/PHYLOGENY

# downloaded reference fastas from ENA/ncbi
KT449822.1_Trichuris_suis.fa
KT449823.1_Trichuris_suis.fa
KT449824.1_Trichuris_sp._baboon.fa
KT449825.1_Trichuris_sp.TTB2.fa
KT449826.1_Trichuris_trichiura.fa
NC_002681.1_Trichinella_spiralis.fa
NC_017747.1_Trichuris_suis.fa
NC_017750.1_Trichuris_trichiura.fa
NC_018596.1_Trichuris_discolor.fa
NC_018597.1_Trichuris_ovis.fa
NC_028621.1_Trichuris_muris.fa


ln -s ../../04_VARIANTS/GATK_HC_MERGED/mito_samples3x_missing0.8.recode.vcf

bgzip mito_samples3x_missing0.8.recode.vcf
tabix -p vcf mito_samples3x_missing0.8.recode.vcf.gz

bcftools query --list-samples mito_samples3x_missing0.8.recode.vcf > samples.list

REF=../../01_REF/trichuris_trichiura.fa
MITO_REF=Trichuris_trichiura_MITO
VCF=mito_samples3x_missing0.8.recode.vcf.gz

while read SAMPLE; do \
     samtools faidx ${REF} ${MITO_REF} | bcftools consensus ${VCF} --missing N --output ${SAMPLE}.mito.fa --sample ${SAMPLE};
done < samples.list

# # mask to only include protein coding regions
# while read SAMPLE; do \
#      samtools faidx ${REF} ${MITO_REF} | bcftools consensus ${VCF} --mask /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/01_REF/mito.mask.bed --missing N --output ${SAMPLE}.mito.fa --sample ${SAMPLE};
# done < samples.list


# fix sample names in the file
for i in *.mito.fa; do \
     sed -i "s/Trichuris_trichiura_MITO/${i%.mito.fa}/" ${i};
done

# bring all of the pseudoreferences together with the ENA/NCBI references
cat *mito.fa > all_pseudoreferences.fa


module load mafft/7.407=1-c1

bsub.py --threads 20 10 mafft "mafft --thread 20 --maxiterate 1000 --globalpair all_pseudoreferences.fa \> AN_MN_pseudomtDNA.aln"

bsub.py --threads 20 10 mafft_refs "mafft --thread 20 --maxiterate 1000 --globalpair refs.fa \> refs.aln"

```

```bash
for i in *trimmed.bam; do
     bcftools mpileup -C 50 --min-MQ 20 --min-BQ 30 --skip-indels ${i} --fasta-ref ../01_REF/trichuris_trichiura.fa --regions Trichuris_trichiura_MITO | bcftools call --ploidy 1 -c -Oz -o ${i%.trimmed.bam}.vcf.gz

     tabix ${i%.trimmed.bam}.vcf.gz

     bedtools complement -i ${i%.trimmed.bam}.vcf.gz -g mito.genome > mask.bed

     bcftools consensus --mask mask.bed --missing N -f ../01_REF/Trichuris_trichiura_MITO.fasta ${i%.trimmed.bam}.vcf.gz -o ${i%.trimmed.bam}.fa

     sed -i "s/Trichuris_trichiura_MITO/${i%.trimmed.bam}/g" ${i%.trimmed.bam}.fa;
done

```
