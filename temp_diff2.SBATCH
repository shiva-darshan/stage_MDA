#!/bin/bash
# File temp_diff2.SBATCH
#SBATCH --nodes=1
#SBATCH -c 9
#SBATCH --job-name="FPUT_chain_therm_tension"
#SBATCH --output=%x%J.out
#SBATCH --mail-user=shiva.darshan@universite-paris-saclay.fr
#SBATCH --mail-type=BEGIN,FAIL,END

echo Host `hostname`

### Nombre total de CPU
echo Il a ete aloue $SLURM_JOB_CPUS_PER_NODE cpus

/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 0 0 100 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 0.25 0 100 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 0.5 0 100 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 0.75 0 100 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 1 0 100 0.005 2000000&

/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 0 0 100 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 0.25 0 100 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 0.5 0 100 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 0.75 0 100 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 1 0 100 0.005 2000000; echo finished&