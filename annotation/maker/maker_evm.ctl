#-----Transcript weights
evmtrans=3 #default weight for source unspecified est/alt_est alignments
evmtrans:blastn=3 #weight for blastn sourced alignments
evmtrans:est2genome=10 #weight for est2genome sourced alignments
evmtrans:tblastx=3 #weight for tblastx sourced alignments
evmtrans:cdna2genome=10 #weight for cdna2genome sourced alignments

#-----Protein weights
evmprot=5 #default weight for source unspecified protein alignments
evmprot:blastx=5 #weight for blastx sourced alignments
evmprot:protein2genome=10 #weight for protein2genome sourced alignments
evmprot:gth=10

#-----Abinitio Prediction weights
evmab=1 #default weight for source unspecified ab initio predictions
evmab:snap=2 #weight for snap sourced predictions
evmab:augustus=10 #weight for augustus sourced predictions
evmab:fgenesh=10 #weight for fgenesh sourced predictions
evmab:GeneMark.hmm=5
evmab:Genscan=5
evmab:Genewise=10
evmab:GlimmerHMM=2


