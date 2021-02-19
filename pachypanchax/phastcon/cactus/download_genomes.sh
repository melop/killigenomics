
lftp -e 'mirror --parallel=4 -c -L ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/266/775/GCA_001266775.1_Austrofundulus_limnaeus-1.0 austrofundulus; exit' -u rcui@age.mpg.de, ftp.ncbi.nlm.nih.gov
lftp -e 'mirror --parallel=4 -c -L ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/775/205/GCA_002775205.2_X_maculatus-5.0-male xmac5.0male; exit' -u rcui@age.mpg.de, ftp.ncbi.nlm.nih.gov
lftp -e 'mirror --parallel=4 -c -L ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/826/765/GCA_000826765.1_Fundulus_heteroclitus-3.0.2 Fhet3.0.2; exit' -u rcui@age.mpg.de, ftp.ncbi.nlm.nih.gov
lftp -e 'mirror --parallel=4 -c -L ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/649/575/GCA_001649575.1_ASM164957v1 kryptolebias; exit' -u rcui@age.mpg.de, ftp.ncbi.nlm.nih.gov



