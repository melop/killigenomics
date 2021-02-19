source /beegfs/group_dv/software/source/progressiveCactus2/progressiveCactus/environment
mkdir split_NFZ_maf
cd split_NFZ_maf
hal2mafMP.py --splitBySequence --smallSize 10000000 --numProc 30 --noAncestors --noDupes --refGenome NFZ --targetGenomes NFZ,NOR,AAU,CTO,PLP,kryptolebias,austrofundulus ../progressiveCactusAln_9spp.hal NFZ_ref_5spp.maf

