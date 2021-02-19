cat /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/exonerate/NFZ_v2.0/pipeline/sum_longest_transcript.gff | awk '{if ($3=="CDS") print "NFZ."$0 }' > NFZ2.0/cds.gff
