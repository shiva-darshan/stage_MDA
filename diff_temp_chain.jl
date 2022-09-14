using JLD

include("flux.jl")
include("atom_chains.jl")


function main(pot_type, a, b, alpha, d, T_l, T_r, tension, lambda, T, burn_in_time, time_step, 
	save_cov_mat)
	gamma_l, gamma_r = 2*lambda, 2*lambda
	T_mid = (T_l + T_r)/2
	F(t) = 0

	left_atom_fixed = false

	if pot_type == "FPUT"
		U, grad_U, v, v_prime = get_pot_fput(alpha, a, b)
		#for now rejection sampling only works for a,b = -2.5, 1 and T >= 0.33
		q_0 = reject_sampling_fput(d, v, T_mid) 
		p_0 = sqrt(T_mid) * randn(d)

		rotar = false
		recenter = alpha == 0

	elseif pot_type == "harmonic"
		U, grad_U, v, v_prime = get_pot_harm(alpha, a)
		q_0 = sqrt(T_mid) * randn(d)
		p_0 = sqrt(T_mid) * randn(d)

		rotar = false
		recenter = alpha == 0

	elseif pot_type == "cos"
		U, grad_U, v, v_prime = get_pot_cos(alpha, a)
		q_0 = reject_sampling_cos(d, v, T_mid)
		p_0 = sqrt(T_mid) * randn(d)

		rotar = true
		recenter = false

	else
		throw(ArgumentError(string("\"", pot_type, "\" is not a valid potential type argument", 
            "\nValid potential type arguments: \"FPUT\", \"harmonic\", \"cos\"")))
	end


	N = round(Int, T/time_step)

	if save_cov_mat
		averages, T_profile, cov_mat, r_profile, e_profile, e_flux, p_flux = atom_chain(q_0, p_0, 
		grad_U, v, v_prime, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, T, N, burn_in_time, 
		save_cov_mat; rotar = rotar, pinning = alpha, recenter = recenter, left_atom_fixed = left_atom_fixed)

		save(string("./data/", pot_type, "_fluxes_a", a, "_b", b, "_alpha", alpha, "_Tl", T_l, 
			"_Tr", T_r, "_lambda", lambda, "_tension", tension, "_d", d, "_T", T, "_burn", 
			burn_in_time, "_dt", time_step, ".jld"), "p_flux", p_flux, "e_flux", e_flux)

		save(string("./data/", pot_type, "_profiles_w_cov_a", a, "_b", b, "_alpha", alpha, "_Tl", T_l, 
			"_Tr", T_r, "_lambda", lambda, "_tension", tension, "_d", d, "_T", T, "_burn", 
			burn_in_time, "_dt", time_step, ".jld"), "T_profile", T_profile, "cov_mat", cov_mat, 
			"r_profile", r_profile, "e_profile", e_profile, "averages", averages)

	else
		averages, T_profile, q_var, r_profile, e_profile, e_flux, p_flux = atom_chain(q_0, p_0, 
		grad_U, v, v_prime, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, T, N, burn_in_time; 
		rotar = rotar, pinning = alpha, recenter = recenter, left_atom_fixed = left_atom_fixed)

		save(string("./data/", pot_type, "_fluxes_a", a, "_b", b, "_alpha", alpha, "_Tl", T_l, 
			"_Tr", T_r, "_lambda", lambda, "_tension", tension, "_d", d, "_T", T, "_burn", 
			burn_in_time, "_dt", time_step, ".jld"), "p_flux", p_flux, "e_flux", e_flux)

		save(string("./data/", pot_type, "_profiles_a", a, "_b", b, "_alpha", alpha, "_Tl", T_l, 
			"_Tr", T_r, "_lambda", lambda, "_tension", tension, "_d", d, "_T", T, "_burn", 
			burn_in_time, "_dt", time_step, ".jld"), "T_profile", T_profile, "q_var", q_var, 
			"r_profile", r_profile, "e_profile", e_profile, "averages", averages)
	end
end


if abspath(PROGRAM_FILE) == @__FILE__
	pot_type = ARGS[1]
	a = parse(Float64, ARGS[2])
	b = parse(Float64, ARGS[3])
	alpha = parse(Float64, ARGS[4])
	T_l = parse(Float64, ARGS[5])
	T_r = parse(Float64, ARGS[6])
	tension = parse(Float64, ARGS[7])
	lambda = parse(Float64, ARGS[8])
	d = parse(Int64, ARGS[9])
	dt = parse(Float64, ARGS[10])
	T = parse(Float64, ARGS[11])
	burn_in_time = parse(Float64, ARGS[12])
	save_cov_mat = parse(Bool, ARGS[13])

	print(string("Started: Potential = ", pot_type, ", a = ", a, ", b = ", b, ", alpha = ", 
				alpha, ",\nTl = ", T_l, ", Tr = ", T_r, ", lambda = ", lambda, ", tension = ", 
				tension, ",\nd = ", d, ", T = ", T, ", burn-in time = ", burn_in_time, ", time step = ", 
				dt, ",\nSaving Covariance Matrix = ", save_cov_mat, "\n\n"))
	flush(stdout)
	#some stuff still relies on the hardcoded parameters for the potential
	main(pot_type, a, b, alpha, d, T_l, T_r, tension, lambda, T, burn_in_time, dt, save_cov_mat)
	print(string("Finished: Potential = ", pot_type, ", a = ", a, ", b = ", b, ", alpha = ", 
				alpha, ",\nTl = ", T_l, ", Tr = ", T_r, ", lambda = ", lambda, ", tension = ", 
				tension, ",\nd = ", d, ", T = ", T, ", burn-in time = ", burn_in_time, ", time step = ", 
				dt, ",\nSaving Covariance Matrix = ", save_cov_mat, "\n\n"))
	flush(stdout)
end