PIDS=()
a=1
b=0
alpha=0
T_l=0.3
gamma_l=2
tension=0
lambda=0
d=1000
dt=0.005
T=100000000
burn=10000
period=0
#omega=1
save_cov=false
c=20
omegas=(1 2 3 4 5 6 7 8 9 10)

for omega in ${omegas[@]}; do
	#pot_type, a, b, c, alpha, omega, T_l, gamma_l, tension, lambda, d, time_step, T, burn_in_time, recording_period, save_cov
	/home/users/darshan/julia-1.7.2/bin/julia forced_chain.jl cos $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $period $save_cov&
	PIDS+=($!)
done

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished