#!/bin/bash
# File temp_diff.SBATCH
#SBATCH --nodes=1
#SBATCH -c 10
#SBATCH --job-name="TEMP_DIFF_CHAIN"
#SBATCH --output=%x%J.out
#SBATCH --mail-user=shiva.darshan@universite-paris-saclay.fr
#SBATCH --mail-type=BEGIN,FAIL,END

echo Host `hostname`

### Nombre total de CPU
echo Il a ete aloue $SLURM_JOB_CPUS_PER_NODE cpus


/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 100 0.005 1000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 200 0.005 2000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 300 0.005 3000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 400 0.005 4000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 500 0.005 5000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 600 0.005 6000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 700 0.005 7000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 800 0.005 8000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 900 0.005 9000000&
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain2.jl harmonic 1 0 0 0.8 1.2 0 1 1000 0.005 10000000; echo finished&

