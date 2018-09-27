Map reads to the reference genome.

1. Run ln.sh to prepare the reference genome for mapping.
2. Run bash map.local.sbatch to map low-coverage reads, skipping deduplication because complexity is high at low coverage
3. Run bash map_hicov.local.sbatch with deduplication on highly covered samples.
4. Proceed to allelefreqxxxxx/ to estimate allele frequencies
