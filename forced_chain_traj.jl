using JLD

include("atom_chains.jl")


function main(pot_type, a, b, c, alpha, omega, d, T_l, gamma_l, lambda, T, burn_in_time, time_step)
	gamma_r = 0.
	T_r = 0.

	F = get_forcing_cos(c, omega)

	if pot_type == "FPUT"
		U, grad_U, v, v_prime = get_pot_fput(alpha, a, b)
		#for now rejection sampling only works for a,b = -2.5, 1 and T >= 0.33
		q_0 = reject_sampling_fput(d, v, T_l) 
		p_0 = sqrt(T_l) * randn(d)

		rotar = false
		recenter = alpha == 0

	elseif pot_type == "harmonic"
		U, grad_U, v, v_prime = get_pot_harm(alpha, a)
		q_0 = sqrt(T_l) * randn(d)
		p_0 = sqrt(T_l) * randn(d)

		rotar = false
		recenter = alpha == 0

	elseif pot_type == "cos"
		U, grad_U, v, v_prime = get_pot_cos(alpha, a)
		q_0 = reject_sampling_cos(d, v, T_l)
		p_0 = sqrt(T_l) * randn(d)

		rotar = true
		recenter = false

	else
		throw(ArgumentError(string("\"", pot_type, "\" not a valid potential type argument", 
            "\nValid potential type arguments: \"FPUT\", \"harmonic\", \"cos\"")))
	end


	N = floor(Int, T/time_step)

	traj = atom_chain_traj(q_0, p_0, grad_U, F, T_l, T_r, gamma_l, gamma_r, lambda, T, N, 
		burn_in_time; rotar = rotar, recenter = recenter)

	save(string("./data/", pot_type, "_traj", "_a", a, "_b", b, "_c", c, 
		"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
		"_d", d, "_T", T, "_dt", time_step, ".jld"), "trajs", 
	traj)

end


if abspath(PROGRAM_FILE) == @__FILE__

	pot_type = ARGS[1]
	a = parse(Float64, ARGS[2])
	b = parse(Float64, ARGS[3])
	c = parse(Float64, ARGS[4])
	alpha = parse(Float64, ARGS[5])
	omega = parse(Float64, ARGS[6])
	T_l = parse(Float64, ARGS[7])
	gamma_l = parse(Float64, ARGS[8])
	lambda = parse(Float64, ARGS[9])
	d = parse(Int64, ARGS[10])
	dt = parse(Float64, ARGS[11])
	T = parse(Float64, ARGS[12])
	burn_in_time = parse(Float64, ARGS[13])


	print(string("Started: Trajs, Potential = ", pot_type, ", a = ", a, ", b = ", b, ", c = ", c, 
		", alpha = ", alpha, ", omega = ", omega, ",\n",
		"T_l = ", T_l, ", gamma_l = ", gamma_l, ", lambda = ", lambda, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, "\n\n"))
	flush(stdout)


	main(pot_type, a, b, c, alpha, omega, d, T_l, gamma_l, lambda, T, burn_in_time, dt)
	print(string("Finished: Trajs Potential = ", pot_type, ", a = ", a, ", b = ", b, ", c = ", c, 
		", alpha = ", alpha, ", omega = ", omega, ",\n",
		"T_l = ", T_l, ", gamma_l = ", gamma_l, ", lambda = ", lambda, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, "\n\n"))
	flush(stdout)
end


