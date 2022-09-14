PIDS=()
alpha=0
T_l=0.1
gamma_l=2
gamma_r=2
tension=0
lambda=1
dt=0.005
burn=1000000
omega=1
skip=200000
d=100
T=250000000
seed=9080822348
c=10


# global min at x = 1
a=-5
b=4
T_r=4.058
/home/users/darshan/julia-1.7.2/bin/julia forced_transport_coef.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia thermalized_transport_coef.jl FPUT $a $b $alpha $T_l $T_r $gamma_l $gamma_r $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)

# global min at x = 2
a=-2.5
b=1
T_r=5.625
/home/users/darshan/julia-1.7.2/bin/julia forced_transport_coef.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia thermalized_transport_coef.jl FPUT $a $b $alpha $T_l $T_r $gamma_l $gamma_r $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)

# global min at x = 3
a=-1.6667
b=0.4444
T_r=6.807
/home/users/darshan/julia-1.7.2/bin/julia forced_transport_coef.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia thermalized_transport_coef.jl FPUT $a $b $alpha $T_l $T_r $gamma_l $gamma_r $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)


# global min at x = -1
a=5
b=4
T_r=4.057
/home/users/darshan/julia-1.7.2/bin/julia forced_transport_coef.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia thermalized_transport_coef.jl FPUT $a $b $alpha $T_l $T_r $gamma_l $gamma_r $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)


# global min at x = -2
a=2.5
b=1
T_r=5.6
/home/users/darshan/julia-1.7.2/bin/julia forced_transport_coef.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia thermalized_transport_coef.jl FPUT $a $b $alpha $T_l $T_r $gamma_l $gamma_r $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)


# global min at x = -3
a=1.6667
b=0.4444
T_r=6.797
/home/users/darshan/julia-1.7.2/bin/julia forced_transport_coef.jl FPUT $a $b $c $alpha $omega $T_l $gamma_l $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia thermalized_transport_coef.jl FPUT $a $b $alpha $T_l $T_r $gamma_l $gamma_r $tension $lambda $d $dt $T $burn $skip $seed&
PIDS+=($!)



for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished