PIDS=()
PIDS1=()
PIDS2=()

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
T=25000000
save_cov=false
seed=9080822348

c=10

as=(-5 -2.5 -1.6667 1.6667 2.5 5)
bs=(4 1 0.4444 0.4444 1 4)


for ((j=0; j < ${#bs[@]}; j++)); do
	julia forced_chain.jl FPUT ${as[$j]} ${bs[$j]} $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $period $save_cov $seed&
	PIDS+=($!)
done

new_T=25000000
new_burn=0
cont=0

for ((j=0; j < ${#bs[@]}; j++)); do 
	wait ${PIDS[$j]} 
	echo done ${PIDS[$j]}
	julia cont_forced_chain.jl FPUT ${as[$j]} ${bs[$j]} $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $period $save_cov $seed $new_T $new_burn $cont&
	PIDS1+=($!)
done
echo finished run 0


new_T2=50000000
new_burn2=0
cont=1

for ((j=0; j < ${#bs[@]}; j++)); do 
	wait ${PIDS1[$j]} 
	echo done ${PIDS1[$j]}
	julia cont_forced_chain.jl FPUT ${as[$j]} ${bs[$j]} $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $new_T $new_burn $period $save_cov $seed $new_T2 $new_burn2 $cont&
	PIDS2+=($!)
done
echo finished run 1

for i in ${PIDS2[@]}; do 
	wait $i 
	echo done $i
done
echo finished