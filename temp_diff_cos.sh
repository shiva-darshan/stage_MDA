PIDS=()
a=1
b=0
alpha=0
T_l=1
tension=0
lambda=1
d=100
dt=0.005
T=20000000
burn=10000
save_cov=false

T_rs=(1.0877 1.1211 1.1599 1.2397 1.3292 1.4119 1.5339 1.6524 1.8559 2.0699 2.3517 2.7121 3.0697 3.7158 4.6521 5.8193)

for T_r in ${T_rs[@]}; do
	#pot_type, a, b, alpha, T_l, T_r, tension, lambda, d, dt, T, burn_in_time, save_cov
	/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl cos $a $b $alpha $T_l $T_r $tension $lambda $d $dt $T $burn $save_cov&
	PIDS+=($!)
done

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished