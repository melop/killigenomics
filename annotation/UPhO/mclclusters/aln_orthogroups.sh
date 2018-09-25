mkdir UPhO_Seqs
Get_fasta_from_Ref.py -o UPhO_Seqs -q ClustR_m2.txt -R ../allproteins.fa
cd UPhO_Seqs/
bash ../paMATRAX_ray.sh -s -c 
