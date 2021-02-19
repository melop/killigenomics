#!/usr/bin/bash
#SBATCH  -p long
#SBATCH -c 40


/beegfs/group_dv/software/source/G-PhoCS/bin/G-PhoCS PLP.ctl -n 40
