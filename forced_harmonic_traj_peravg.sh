PIDS=()
a=1
b=0
alpha=1
T_l=1
gamma_l=1
lambda=1
d=20
dt=0.005
T=2000000
burn=10000
period=100
sampling=20
per_inds=all
traj_inds=5,10,15,20

cs=(8 16)
omegas=(1)

for c in ${cs[@]}; do
	for omega in ${omegas[@]}; do
		/home/users/darshan/julia-1.7.2/bin/julia forced_chain_traj.jl harmonic $a $b $c $alpha $omega $T_l $gamma_l $lambda $d $dt $T $burn $sampling $period $per_inds $traj_inds&
		PIDS+=($!)
		
	done
done

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished