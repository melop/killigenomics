#converting format
#python2 /beegfs/group_dv/software/source/genomics_general/VCF_processing/parseVCF.py --minQual 60 --skipIndels -i PLP_merged.filtered.vcf.gz -o PLP_merged.filter.geno.gz
python2 /beegfs/group_dv/software/source/genomics_general/popgenWindows.py -g PLP_merged.filter.geno.gz -o PLP_merged.Fst.Dxy.pi.filter.csv.gz \
   -f phased -w 20000 -m 5000 -s 20000 \
   -p PLPA -p PLPB -p PLPC -p PLPD -p PLPE -p PLPF -p PLPG \
   --popsFile pop_file.txt \
   -T 20


