Fix maker gene annotations
1. autocorrect.sh conducts a few steps of fixing automatically. Edit this file and configure the paths.
2. The last step of autocorrect.sh produces maker.replace_frag_aln.gff and all fasta files for the previous steps.
3. Pick up from here manually, run removeoverlap_final.sh to produce maker.replace_frag_aln.noolverlap.gff.
4. Run addgenesymbol.php to assign gene symbols, producing maker.genesymbol.gff
5. Run addgenesymbol.supplement.php to refine gene symbols, producing maker.genesymbol.supple.gff
6. Run finalannotation.php, then run finalannot_sum.php to summarize the runs, producing maker.finalannot.gff
7. Run improve_frag_models.sbatch, followed by improve_frag_models.sum.php to produce maker.finalannot.improved.gff
8. Use commands in maker.refine.proteins.sh to export proteins and mRNA sequences
9. whitelistgenes.php outputs a white list for CAFE analysis
10. plotqual.sh contains plotting commands to visualize fixing results by calling the exportquality.php command and plotqual.R script

