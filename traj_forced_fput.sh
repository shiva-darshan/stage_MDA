PIDS=()
a=-2.5
b=1
alpha=1
T_l=1
gamma_l=2
lambda=1
d=100
dt=0.005
T=20000
burn=0
period=100
omega=1
save_cov=false

cs=(0.5)

for c in ${cs[@]}; do
	#pot_type, a, b, c, alpha, omega, T_l, gamma_l, lambda, d, time_step, T, burn_in_time, recording_period, save_cov
	julia forced_chain_traj.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $lambda $d $dt $T $burn&
	PIDS+=($!)
done

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished