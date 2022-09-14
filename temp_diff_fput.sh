PIDS=()
alpha=0
T_l=0.1
gamma_l=2
tension=0
lambda=1
dt=0.005
burn=1000000
period=0
omega=1
d=100
T=250000000
save_cov=false
seed=9080822348



julia diff_temp_chain.jl FPUT $a $b $alpha $T_l $T_r $tension $lambda $d $dt $T $burn $save_cov&
PIDS+=($!)


for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished