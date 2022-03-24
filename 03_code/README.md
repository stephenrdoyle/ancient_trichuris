# Population genomics of ancient and modern *Trichuris trichiura*

## Overview of workflow and code used in the analysis

### [Genome analyses](ancient_trichuris.00_genome_analyses.md)
- Reference genome
- Genome completeness using BUSCO
- Transfer of gene models using LiftOff
- Running interproscan on liftoff annotation

### [Mapping](ancient_trichuris.01_mapping.md)
- Project setup
- Prepare reference genome
- Raw sequence data
- metadata
- trimming of the raw reads
- FIX: merge duplicate ancient read sets
- Mapping

### [Kraken](ancient_trichuris.02_kraken_contmaination.md)
- kraken - modern samples
- kraken - ancient samples

### [DNA damage](ancient_trichuris.03_dna_damage.md)
- pmdtools
- plot deamination frequencies
- trim bases from reads in bam to remove deamination

### [Genome coverage](ancient_trichuris.04_genome_coverage.md)
- Genome wide coverage
- Generate quantitative stats on coverage for supplementary tables
- generate some coverage plots
- Genome-wide coverage to determine worm sex

### [Variant calling and filtering](ancient_trichuris.05_variant_calling_and_filtering.md)
- Genome scope to estimate heterozyosity
- GATK
- Filter the VCF - SNPable
- Filter the VCF - hardfilter
- Querying SNP and INDEL QC profiles to determine thresholds for filters
- Applying filters to the variants
- merge VCFs
- Filter genotypes based on depth per genotype
- Sample missingness
- Generate an ALL SITES variant set for running pixy

### [Sampling sites](ancient_trichuris.06_sampling_site_maps_and_data.md)
- World map
- Sampling timepoints
- Map of ancient DNA sites, with Fst data showing degree of connectivity between populations

### [PCAs](ancient_trichuris.07_PCAs.md)
- PCA of mitochondrial variants
- PCA of nuclear variants

### [ANGSD](ancient_trichuris.08_ANGSD.md)
- IBS and Coviance matrices
- Single allele analyses of mitochondrial and nuclear variants

### [NGSadmix](ancient_trichuris.09_NGSadmix.md)
- Admixture plots
- Clumpak to determine optimal K

### [Treemix](ancient_trichuris.10_treemix.md)
- Treemix
- plotting treemix data
- Estimating the optimal number of migration edges
- calculate the variance explained by the data

### [Admixtools](ancient_trichuris.11_admixtools.md)
- Prepared data and run admixtools
- plotting admix data

### [Population history using SMC++](ancient_trichuris.12_smc++.md)
- run smc++
- plot all smcpp plots together

### [Genome-wide genetic variation](ancient_trichuris.13_genomewide_genetic_variation.md)
- Running pixy to calculate nucleotide diversity, dXY and Fst between groups
- analyses of nucleotide diversity (Pi)
- dXY and Fst
- extracting top X% of Fst values for each comparison
- Private and shared variation between populations
- Ancient DNA - nucleotide diversity vs sampling age
- Sample Heterozygosity vs sequencing coverage
- Relatedness and kinship between samples in a population
- genome-wide analyses presented in the original draft of the manuscript

### [Beta-tubulin variation](ancient_trichuris.14_beta-tubulin.md)
- Analysis of variation within beta tubulin
- Reanalysis of nucleotide diversity around beta-tubulin


### Original sampling data
- Modern genomes  from human hosts (49)
     - Human worm isolates (27)
          - Uganda (12), China (7), Honduras (8)
     - Human egg isolates (10)
          - Cameroon (5), Tanzania (5)

- Modern genomes from non-human hosts (7)
     - Worm isolates (7)
          - Baboon – Denmark (2), Colobus - Spain (2), Leaf-Monkey – China (3)

- Ancient samples (17 samples from 10 sites)
     - 1000 AD: Viborg, DK (2) – 2x Mito + 0,01x/0,16x NG
     - 1350 AD : Kampen, Netherlands (2) – 2x Mito + 1,8x/4,7x NG - yes
     - 1350-1400 AD: Odense, DK (4) – in analysis - yes
     - 1500-1600 AD: Zwolle, Netherlands (1) – Mito + 2,3x NG - yes
     - 1550-1580 AD: Vilnius, Lithuania (1) – in analysis - yes
     - 1600-1800 AD: Gl. Strand, CPH, DK (4) – in analysis - yes
     - 1680 AD: Kultorvet, CPH, DK (2) – in analysis - yes
     - 1700 AD: Adelgade, CPH, DK (1) – in analysis - yes
