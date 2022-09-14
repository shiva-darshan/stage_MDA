PIDS=()
a=1
b=0
alpha=1
T_l=1
gamma_l=2
lambda=1
d=100
dt=0.005
T=5000000
burn=10000
period=1
omega=1
q_upper=5
q_lower=-5
q_dx=0.01
p_upper=5
p_lower=-5
p_dx=0.01
r_upper=5
r_lower=-5
r_dx=0.01

cs=(0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)

for c in ${cs[@]}; do
	#pot_type, a, b, c, alpha, omega, T_l, gamma_l, lambda, d, time_step, T, burn_in_time, recording_period, period, q_upper, q_lower, q_dx, p_upper, p_lower, p_dx, r_upper, r_lower, r_dx
	/home/users/darshan/julia-1.7.2/bin/julia forced_chain_hist.jl harmonic $a $b $c $alpha $omega $T_l $gamma_l $lambda $d $dt $T $burn $period $q_upper $q_lower $q_dx $p_upper $p_lower $p_dx $r_upper $r_lower $r_dx&
	PIDS+=($!)
done

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished