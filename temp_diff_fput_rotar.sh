PIDS=()
#pot_type, a, b, alpha, T_l, T_r, tension, lambda, d, dt, T, burn_in_time, sampling_interval, inds_to_save
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 0 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 0.25 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 0.5 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 0.75 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.8 1.2 1 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 0 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 0.25 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 0.5 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 0.75 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl FPUT -2.5 1 0 0.9 1.1 1 1 100 0.005 2000000 10000 1 false&
PIDS+=($!)

/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl cos 1 0 0 0.8 1.2 0.75 0 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl cos 1 0 0 0.9 1.1 0.25 0 100 0.005 2000000 10000 1 false&
PIDS+=($!)
/home/users/darshan/julia-1.7.2/bin/julia diff_temp_chain.jl cos 1 0 0 0.9 1.1 0.5 0 100 0.005 2000000 10000 1 false&
PIDS+=($!)

for i in ${PIDS[@]}; do 
	wait $i
	echo done $i
done
echo finished