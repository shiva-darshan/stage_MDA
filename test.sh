# a=1
# b=0
alpha=0
tension=1.
T_l=1
gamma_l=2
lambda=1
d=100
dt=0.005
T=200
burn=10
period=0
omega=1
save_cov=false
seed=1
c=10

as=(-5 -2.5 -1.6667 1.6667 2.5 5)
bs=(4 1 0.4444 0.4444 1 4)


for ((j=0; j < ${#bs[@]}; j++)); do
	julia forced_chain.jl FPUT ${as[$j]} ${bs[$j]} $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $period $save_cov $seed&
done