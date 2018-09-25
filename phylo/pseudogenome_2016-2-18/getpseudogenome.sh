# Getting alignments using genome reference for concatenated analyses and bucky (allowing missing data):
# -w 1000 break off threshold 1000bp; if more than 1000bp gap exists between two alignment blocks, break them into two genes
# -c list containing references to ref genome and the bcf files
# -R character requirements
# -s 0.5: allow 50% missing per column
# -l 500, final alignment length at least 500bp
# -d 5, coverage < 5x turned to "N"
# -N no,  when applying the final divergence filter, don't count "N" as a difference.
# -o output folder

# -D no = turn off divergence filter. -n include N positions in the ref sequence. -I yes: code indels or multiple non-ref alleles as N, so that coordinate identical with ref.
# -r ratio of reads supporting the alternative allele, for example, 1 means 0 reads supporting ref, 1 supporting alt, or 0 to 2, 0 to 3 etc. if doesn't pass this, masked as N

mkdir ./testaln
hhvm getAlignmentAllowN.php -n yes -q 1 -D no -w 1000000000 -c alltaxa.txt -R charRequirements.txt -s 1 -l 10 -d 1 -N no -r 1 -o ./testaln -I yes -b no
[ $? -ne 0 ] &&  exit 1 ;



arrInput=( `cut -f1 alltaxa.txt` )



#mkdir $outFolder
for i in "${arrInput[@]}" 
do 
echo converting $i ...

hhvm phy2fasta.php -L ${i} -o ${i}_pseudogenome2.fasta &

done
wait
