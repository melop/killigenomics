#!/bin/bash
#SBATCH -c 40
#SBATCH -p blade

#min identity (local) = 0.95, -G 0 local identity, -uL unmatched length <=30% of the longer sequence, set for alternatively spliced isoforms -aS the shorter sequence must be mostly encapsulated (90%) 
cd-hit-est -i Trinity.fasta.transdecoder.mRNA -o Trinity.cdhitout.fasta -c 0.95 -M 200000 -T 40 -G 0 -uL 0.3 -aS 0.9 -d 0


