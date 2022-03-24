# Overview of workflow and code used in the analysis


### Genome analyses
- Reference genome
- Genome completeness using BUSCO
- Transfer of gene models using LiftOff
- Running interproscan on liftoff annotation

### Mapping
- Project setup
- Prepare reference genome
- Raw sequence data
- metadata
- trimming of the raw reads
- FIX: merge duplicate ancient read sets
- Mapping

### Kraken
- kraken - modern samples
- kraken - ancient samples

### DNA damage
- pmdtools
- plot deamination frequencies
- trim bases from reads in bam to remove deamination

### Genome coverage
- Genome wide coverage
- Generate quantitative stats on coverage for supplementary tables
- generate some coverage plots
- Genome-wide coverage to determine worm sex

### Variant calling and filtering
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

### Sampling sites
- World map
- Sampling timepoints
- Map of ancient DNA sites, with Fst data showing degree of connectivity between populations

### PCAs
- PCA of mitochondrial variants
- PCA of nuclear variants

### ANGSD
- IBS and Coviance matrices
- Single allele analyses of mitochondrial and nuclear variants

### NGSadmix
- Admixture plots
- Clumpak to determine optimal K

### Treemix
- Treemix
- plotting treemix data
- Estimating the optimal number of migration edges
- calculate the variance explained by the data

### Admixtools
- Prepared data and run admixtools
- plotting admix data

### Population history using SMC++
- run smc++
- plot all smcpp plots together

### Genome-wide genetic variation
- Running pixy to calculate nucleotide diversity, dXY and Fst between groups
- analyses of nucleotide diversity (Pi)
- dXY and Fst
- extracting top X% of Fst values for each comparison
- Private and shared variation between populations
- Ancient DNA - nucleotide diversity vs sampling age
- Sample Heterozygosity vs sequencing coverage
- Relatedness and kinship between samples in a population
- genome-wide analyses presented in the original draft of the manuscript

### Beta-tubulin variation
- Analysis of variation within beta tubulin
- Reanalysis of nucleotide diversity around beta-tubulin
