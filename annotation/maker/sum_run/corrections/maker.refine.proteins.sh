hhvm gff2protein.php -R query.fa -M split__main.gff -o maker.proteins.fa -t maker.transcripts.fa
hhvm gff2protein.php -R query.fa -M chimericfixed.gff -o maker.split.proteins.fa -t maker.split.transcripts.fa
hhvm gff2protein.php -R query.fa -M maker.refined.gff -o maker.refined.proteins.fa -t maker.refined.transcripts.fa
hhvm gff2protein.php -R query.fa -M maker.refined.removeisoform.gff -o maker.refined.noisoform.proteins.fa -t maker.refined.noisoform.transcripts.fa
hhvm gff2protein.php -R query.fa -M maker.refined.removeisoform.missedaddedback.gff -o maker.addedback.proteins.fa -t maker.addedback.transcripts.fa
hhvm gff2protein.php -R query.fa -M maker.replace_frag_aln.gff -o maker.frag2aln.proteins.fa -t maker.frag2aln.transcripts.fa
hhvm gff2protein.php -R query.fa -M maker.finalannot.improved.gff -o maker.finalannot.improved.proteins.fa -t maker.finalannot.improved.transcripts.fa

#hhvm gff2protein_ext.php -R query.fa -M maker.refined.removeisoform.gff -o maker.refined.noisoform.trimmed.proteins.fa -t maker.refined.noisoform.trimmed.transcripts.fa

