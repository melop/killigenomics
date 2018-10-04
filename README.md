# killigenomics
Scripts for African killifish genomics. The raw data can be found on NCBI SRA xxxxx.
The main programming language used was PHP >5. They can be most efficiently run with HHVM or PHP 7. 

* genome_assembly/ scripts for performing genome assembly and improvements
* annotation/ pipelines for genome annotation and improvements
* trinity/ scripts for performing de novo RNAseq assembly
* pseudogenomes/ scripts for generating pseudogenomes by switching a SNP based on mapped reads. 
* phylo/ build phylogenetic trees
* assemble_mt/ contains scripts for mitochondrial genome assembly, annotation and analyses
* relax_codeml/ contains scripts for alignment generation, gene tree estimation, RELAX and CODEML analyses, amino acid convergence analysis
* WGS_seq/ contains scripts for processing population resequencing data
* cpp/ contains c++ code that needs compilation
* longitudinalrnaseq/ contains scripts to download and reanalyze a previously published longitudinal RNAseq dataset. 

Third-party software packages include: 

* Trimmomatic 0.32
* Preseq 1.0.2
* NextClip 1.3.1
* Discovar de novo 52488
* NCBI BLAST 2.2.31+
* BESST (git version 7dba5dd)
* BESST_RNA (git version 9cc039b)
* MasurCA 3.2.1
* Gapfiller 1.10
* Reapr 1.0.18
* Redundans 0.13a
* SGA preQC 0.10.13
* BWA-MEM 0.7.12-r1039
* SMALT 0.5.3
* STAR 2.4.2a
* BUSCO78 v2beta4
* Metassembler  1.3
* PILON 1.22
* AllMaps 0.6.9
* ProgressiveCactus (git version 7af8f26)
* caper (R package, 1.0.1)
* RAxML 8.2.4
* TranslatorX 1.1
* MAFFT v7.305b
* Gblocks 0.91b
* Consel 1.20
* NOVOPlasty 2.6.3
* MITOS server (http://mitos.bioinf.uni-leipzig.de/index.py)
* Trinity 2.1.1
* Maker 3.0beta
* Augustus 3.2.1
* CDHitEST 4.6.4
* Exonerate 2.2.0
* RepeatModeler  1.73
* RepeatMasker 1.331
* UPhO (git version 4ec1589)
* MCL 14-137
* Fasttree 2.1.9
* HMMER 2
* GeneWise 2.4.1
* EvidenceModeler 1.1.1
* CAFE' 4.0.1
* TrimAl106 v1.4.rev15
* CIWOG/GECA (git version ede81942)
* SGA preQC 0.10.13
* GATK 3.4.46
* Samtools/Bcftools 1.2
* Picard v.1.119
* APE 5.1
* PHYTOOLS 0.6.44
* HyPhy v.2.2.5/v2.3.8
* PAML 4.8
* PCOC 07022018
* Consurf 1.0.6
* Rate4Site 3.2
* Package-GFE (05/29/15)
* scrm 1.7
* fastStructure 1.0
* PhyloBayes4.1
* MSMC 2
* Bcftools 1.6
* LDHat (git version 81596e2)
* Anavar 1.4
* wgsim 0.3.1-r13
* featureCounts 1.5.0
* Deseq 2
* Angsd 0.918
