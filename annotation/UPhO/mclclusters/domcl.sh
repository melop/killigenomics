ln -s ../mcl_1e*.abc ./
mcxload -abc mcl_1e-10.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o seq.mci -write-tab seq.tab

mcl seq.mci -te 30 -I 2 -use-tab seq.tab 


sed 's/\t/,/g' out.seq.mci.I20  > mciI20.csv
python -c 'import BlastResultsCluster as BRC; BRC.redundant("mciI20.csv", 2)'

python  /software/source/UPhO/Get_fasta_from_Ref.py -q ClustR_m5.txt -R ../allproteins.fa -o ClusteRs
