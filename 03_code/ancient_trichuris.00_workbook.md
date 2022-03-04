# Population genomics of ancient and modern *Trichuris trichiura*

## Background

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


### Overarching questions
- describe new genome?
     - basic assembly stats
- global population structure
     - worldmap
     - PCA (nuclear / mtDNA)
     - structure
     - dispersal - MSMC (difference in popn size between human and animals?)
- admixture
     - from ancient to modern
     - from humans to animals
          - treemix / admixture / D-statistics
- genome-wide diversity
     - similar to oncho popgen paper
          - nucleotide diversity (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5512550/; Fig2)
               - human vs animals (nucleotide div, Tajimas D)
               - Fst (human vs animal)
               - Fst (two clades of human Tt)

- evidence of selection around beta-tubulin?
     - might apply to some modern samples


### Downstream analyses
- relatedness
     - PCA https://speciationgenomics.github.io/pca/
          - whole genomes modern
          - mtDNA genomes all
     - faststructure
- admixture
     - treemix
          - tutorial: https://speciationgenomics.github.io/Treemix/
     - admixtools https://github.com/uqrmaie1/admixtools
          - D statistics
- population history
     - MSMC / SMC++ (https://github.com/popgenmethods/smcpp)
- genome-wide patterns of genetic diversity
     - Pi
     - Tajimas D
     - Fst
