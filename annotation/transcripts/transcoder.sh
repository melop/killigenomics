#!/bin/bash
#SBATCH -p blade,himem,hugemem
#SBATCH -n 1
#SBATCH -c 40
TRANSCODERDIR=/software/source/TransDecoder-2.0.1/
PFAMDB=/source/Trinotate-2.0.2/seqdb/Pfam-A.hmm
UNIREFDB=/source/Trinotate-2.0.2/seqdb/uniprot_uniref90.trinotate.pep
THREADS=40

echo user: $USER $UID
echo host: `hostname`
#exit

pwd

source ~/.bashrc
module load BLAST 

$TRANSCODERDIR/TransDecoder.LongOrfs -m 100 -t Trinity.fasta

echo Running HMMSCAN...
hmmscan --cpu  $THREADS --domtblout pfam.domtblout $PFAMDB Trinity.fasta.transdecoder_dir/longest_orfs.pep

echo Running blastp
blastp -query Trinity.fasta.transdecoder_dir/longest_orfs.pep  -db $UNIREFDB  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads $THREADS > blastp.outfmt6

echo Predicting...
$TRANSCODERDIR/TransDecoder.Predict --retain_long_orfs 1000 -t Trinity.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6




