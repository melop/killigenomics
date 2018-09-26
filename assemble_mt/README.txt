De novo assembly of mitochondrial genomes from short reads and annotation.
The scaffolds from the AAU, CTO, NOR and PLP assemblies identified as the mitochondrion are used as probes. see mtref/

1. Edit list.txt, run assemble_mt.sbatch on SLURM to obtain mitochondrial genomes. output folder structure: out/[REFNAME]/[SPECIESCODE]
2. If you see a file named Circularized_assembly_1_??.fasta in the output folder, it suggests that the assembler is able to close the mt genome, otherwise you may need to manually close the genome.
3. Proceed to mitos_annotation for gene annotation and analyses.

