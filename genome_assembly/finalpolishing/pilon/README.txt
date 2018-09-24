1. Run mapreads.sh to map the reads to the assembly
2. BWA marks all mate pair reads as improperly paired, this is a problem for pilon. Use fixmatepairflag.sh to fix the bam files
3. run pilon.sh to fix small problems on the assembly
