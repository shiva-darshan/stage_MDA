#!/bin/bash
# File temp_diff.SBATCH
#SBATCH --nodes=1
#SBATCH -c 13
#SBATCH --job-name="temp_diff_fput_rotar"
#SBATCH --output=%x%J.out
#SBATCH --mail-user=shiva.darshan@universite-paris-saclay.fr
#SBATCH --mail-type=BEGIN,FAIL,END

echo Host `hostname`

### Nombre total de CPU
echo Il a ete aloue $SLURM_JOB_CPUS_PER_NODE cpus

bash temp_diff_fput_rotar.sh; echo finished



