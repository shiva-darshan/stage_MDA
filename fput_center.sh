PIDS=()
alpha=0
T_l=0.1
gamma_l=2
tension=0
lambda=0
dt=0.005
burn=1000000
omega=1
skip=2000
d=1000
T=100000000
seed=9104973537
cont=0
T_new=0
burn_new=0

c=10


h=(2 0.6667 -1 -2 -0.6667)

# global min at x = 1
a=-5
b=4


/home/users/darshan/julia-1.7.2/bin/julia fput_center.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed $cont $T_new $burn_new&
PIDS+=($!)


# global min at x = 2
a=-2.5
b=1


/home/users/darshan/julia-1.7.2/bin/julia fput_center.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed $cont $T_new $burn_new&
PIDS+=($!)


# global min at x = 3

a=-1.6667
b=0.4444


/home/users/darshan/julia-1.7.2/bin/julia fput_center.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed $cont $T_new $burn_new&
PIDS+=($!)


# global min at x = -1
a=5
b=4


/home/users/darshan/julia-1.7.2/bin/julia fput_center.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed $cont $T_new $burn_new&
PIDS+=($!)


# global min at x = -2
a=2.5
b=1


/home/users/darshan/julia-1.7.2/bin/julia fput_center.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed $cont $T_new $burn_new&
PIDS+=($!)


# global min at x = -3

a=1.6667
b=0.4444


/home/users/darshan/julia-1.7.2/bin/julia fput_center.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed $cont $T_new $burn_new&
PIDS+=($!)



for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished