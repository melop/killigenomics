sAccFile=../motifscan_meme/plp.nfzref.acc.nonCDS.gff
sControlFile=annualconsv/annual.consv.nonCDS.gff 
sort -k1,1 -k4,4n ../motifscan_meme/annualconsv/gene.cds_edges.gff > gene.cds_edges.sorted.gff

bedtools closest  -a gene.cds_edges.sorted.gff -b $sAccFile -D a -id -io | awk '{if ($10!=".") print $0 }' > plp.nfzref.upstreamcds.gff
nLnCount=`cat $sAccFile | wc -l`
echo subsampling $nLnCount lines from control...

rm annualconserved.nfzref.upstreamcds.gff
for i in `seq 0 99`; do
shuf -n $nLnCount $sControlFile | sort -k1,1 -k4,4n > annual.consv.nonCDS.subsample.gff
bedtools closest  -a gene.cds_edges.sorted.gff -b annual.consv.nonCDS.subsample.gff -D a -id -io | awk '{if ($10!=".") print $0 }' >> annualconserved.nfzref.upstreamcds.gff
done

Rscript upstream_dist_compare.R > chosen_dist_cutoff.txt
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-3548.134) print $0 }' > plp.nfzref.upstreamcds.3548.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-177827.9) print $0 }' > plp.nfzref.upstreamcds.177827.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-10000) print $0 }' > plp.nfzref.upstreamcds.10kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-20000) print $0 }' > plp.nfzref.upstreamcds.20kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-30000) print $0 }' > plp.nfzref.upstreamcds.30kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-40000) print $0 }' > plp.nfzref.upstreamcds.40kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-50000) print $0 }' > plp.nfzref.upstreamcds.50kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-60000) print $0 }' > plp.nfzref.upstreamcds.60kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-70000) print $0 }' > plp.nfzref.upstreamcds.70kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-80000) print $0 }' > plp.nfzref.upstreamcds.80kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-90000) print $0 }' > plp.nfzref.upstreamcds.90kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-100000) print $0 }' > plp.nfzref.upstreamcds.100kb.gff



cat plp.nfzref.upstreamcds.gff | awk '{if($20>-20000 && $20<=-10000 ) print $0 }' > plp.nfzref.upstreamcds.10-20kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-30000 && $20<=-20000 ) print $0 }' > plp.nfzref.upstreamcds.20-30kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-40000 && $20<=-30000 ) print $0 }' > plp.nfzref.upstreamcds.30-40kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-50000 && $20<=-40000 ) print $0 }' > plp.nfzref.upstreamcds.40-50kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-60000 && $20<=-50000 ) print $0 }' > plp.nfzref.upstreamcds.50-60kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-70000 && $20<=-60000 ) print $0 }' > plp.nfzref.upstreamcds.60-70kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-80000 && $20<=-70000 ) print $0 }' > plp.nfzref.upstreamcds.70-80kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-90000 && $20<=-80000 ) print $0 }' > plp.nfzref.upstreamcds.80-90kb.gff
cat plp.nfzref.upstreamcds.gff | awk '{if($20>-100000 && $20<=-90000 ) print $0 }' > plp.nfzref.upstreamcds.90-100kb.gff
