#!/bin/bash
# File cuda_test.SBATCH
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --gres=gpu:1
#SBATCH --job-name="TEST"
#SBATCH --output=%x%J.out
#SBATCH --mail-user=shiva.darshan@universite-paris-saclay.fr
#SBATCH --mail-type=BEGIN,FAIL,END

echo Host `hostname`

### Nombre total de CPU
echo Il a ete aloue $SLURM_JOB_CPUS_PER_NODE cpus
 
 PIDS=()

# /home/users/darshan/julia-1.7.2/bin/julia -e 'include("/home/users/darshan/testcpucuda.jl"); cpu_test(q_n, p_n, N, n)'&
# PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia -e 'include("/home/users/darshan/testcpucuda.jl"); cuda_test(q_n, p_n, N, n)'&
PIDS+=($!)

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished