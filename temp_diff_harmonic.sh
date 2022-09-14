PIDS=()
a=1
b=0
alpha=0
T_l=1
lambda=1
tension=0
d=100
dt=0.005
T=20000000
burn=10000
save_cov=false


if [[ $alpha -eq 1 ]]; then 
	T_rs=(1.0656 1.0945 1.1286 1.1679 1.2125 1.2624 1.3175 1.3778 1.4434 1.5143 1.5904 1.6717 1.7583 1.8501 1.9472 2.0495)
elif [[ $alpha -eq 0 ]]; then
	T_rs=(1.0239 1.0344 1.0468 1.0612 1.0774 1.0956 1.1156 1.1376 1.1615 1.1873 1.2150 1.2447 1.2762 1.3096 1.3450 1.3823)
fi


for T_r in ${T_rs[@]}; do
	#pot_type, a, b, alpha, T_l, T_r, tension, lambda, d, dt, T, burn_in_time, save_cov
	/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl harmonic $a $b $alpha $T_l $T_r $tension $lambda $d $dt $T $burn $save_cov&
	PIDS+=($!)
done

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished