#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --job-name="HARMONIC_FORCED"
#SBATCH --output=%x%J.out
#SBATCH --mail-user=shiva.darshan@universite-paris-saclay.fr
#SBATCH --mail-type=BEGIN,FAIL,END

echo Host `hostname`

### Nombre total de CPU
echo Il a ete aloue $SLURM_JOB_CPUS_PER_NODE cpus

bash forced_harmonic.sh; echo finished&