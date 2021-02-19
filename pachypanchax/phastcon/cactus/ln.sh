#link killifish assemblies
sSp=AAU
mkdir $sSp
ln -sf /beegfs/group_dv/home/RCui/killifish_genomes/denovo/release/${sSp}/1.0/mergegap/blastmergegap_cov/summed.blastmerged.fa ./${sSp}/${sSp}.fna
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/repeatmasker/${sSp}_v1.0_blastmergegap/scf.fa.out ./${sSp}/${sSp}.repeatmask.out
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/maker/${sSp}_final_1/run_blastmerge_cov_exonerateOn/corrections/maker.finalannot.improved.gff ${sSp}/${sSp}.gff

sSp=CTO
mkdir $sSp
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/denovo/release/CTO/1.0/mergegap/blastmergegap_cov/final_polish/pilon_out/output/cto.pilon.rename.fa ./${sSp}/${sSp}.fna
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/repeatmasker/${sSp}_v1.0_pilon/scf.fa.out ./${sSp}/${sSp}.repeatmask.out
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/maker/${sSp}_final_1/run_pilon_exonerateOn/corrections/maker.finalannot.improved.gff ${sSp}/${sSp}.gff

sSp=NOR
mkdir $sSp
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/denovo/release/NOR/1.0/mergegap/blastmergegap_cov/final_polish/pilon_out/output/nor.pilon.rename.fa ./${sSp}/${sSp}.fna
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/repeatmasker/${sSp}_v1.0_pilon/scf.fa.out ./${sSp}/${sSp}.repeatmask.out
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/maker/${sSp}_final_1/run_pilon_exonerateOn/corrections/maker.finalannot.improved.gff ${sSp}/${sSp}.gff

sSp=PLP
mkdir $sSp
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/denovo/release/PLP/1.0/mergegap/blastmergegap_cov/final_polish/pilon_out/output/plp.pilon.rename.fa ./${sSp}/${sSp}.fna
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/repeatmasker/${sSp}_v1.0_pilon/scf.fa.out ./${sSp}/${sSp}.repeatmask.out
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/maker/${sSp}_final_1/run_pilon_exonerateOn/corrections/maker.finalannot.improved.gff ${sSp}/${sSp}.gff

sSp=NFZ2.0
mkdir $sSp
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/denovo/meta/NFZ/LG_allpath.jena.stanford.supernova.canu/Both_Crosses.fasta ./${sSp}/${sSp}.fna
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/repeatmasker/NFZ_v2.0/scf.fa.out ./${sSp}/${sSp}.repeatmask.out
ln -sf /beegfs_old/group_dv/home/RCui/killifish_genomes/annotation/exonerate/NFZ_v2.0/pipeline/sum_longest_transcript.gff ${sSp}/${sSp}.gff



