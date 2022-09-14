PIDS=()
a=1
b=0
alpha=0
T_l=1
gamma_l=2
tension=0
lambda=1
d=100
dt=0.005
T=20000000
burn=10000
period=1
omega=1
save_cov=false

cs=(0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)

for c in ${cs[@]}; do
	#pot_type, a, b, c, alpha, omega, T_l, gamma_l, tension, lambda, d, time_step, T, burn_in_time, recording_period, save_cov
	/home/users/darshan/julia-1.7.2/bin/julia forced_chain.jl harmonic $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $period $save_cov&
	PIDS+=($!)
done

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished