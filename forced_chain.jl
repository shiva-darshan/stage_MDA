using JLD

include("flux.jl")
include("atom_chains.jl")


function main(pot_type, a, b, c, alpha, omega, d, T_l, gamma_l, tension, lambda, T, burn_in_time, 
	time_step, recording_period, save_cov_mat, seed)
	gamma_r = 0.
	T_r = 0.

	F = get_forcing_cos(c, omega, d)

	left_atom_fixed = false

	Random.seed!(seed)
	RNG = copy(Random.default_rng())

	if pot_type == "FPUT"
		U, grad_U, v, v_prime, v_second = get_pot_fput(alpha, a, b)
		
		q_0 = reject_sampling_fput(RNG, d, v, a, b, T_l) 
		p_0 = sqrt(T_l) * randn(RNG, d)

		rotar = false
		recenter = alpha==0

	elseif pot_type == "harmonic"
		U, grad_U, v, v_prime, v_second = get_pot_harm(alpha, a)
		q_0 = sqrt(T_l) * randn(RNG, d)
		p_0 = sqrt(T_l) * randn(RNG, d)

		rotar = false
		recenter = alpha==0

	elseif pot_type == "cos"
		U, grad_U, v, v_prime, v_second = get_pot_cos(alpha, a)
		q_0 = reject_sampling_cos(RNG, d, v, T_l)
		p_0 = sqrt(T_l) * randn(RNG, d)

		rotar = true
		recenter = false

	else
		throw(ArgumentError(string("\"", pot_type, "\" not a valid potential type argument", 
            "\nValid potential type arguments: \"FPUT\", \"harmonic\", \"cos\"")))
	end


	N = floor(Int, T/time_step)
	params = Dict([("pot_type", pot_type), ("a", a), ("b", b), ("c", c), ("alpha", alpha), 
		("omega", omega), ("d", d), ("T_l", T_l), ("gamma_l", gamma_l), ("tension", tension), 
		("lambda", lambda), ("T", T), ("burn_in_time", burn_in_time), ("dt", dt), 
		("recording_period", recording_period), ("save_cov_mat", save_cov_mat), ("seed", seed)])

	if save_cov_mat
		if recording_period > 0
			state, period_avg, T_profile, cov_mat, r_profile, e_profile, e_flux, p_flux = atom_chain(q_0, p_0, 
				grad_U, v, v_prime, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, T, N, burn_in_time, 
				recording_period, save_cov_mat; rotar = rotar, pinning = alpha, recenter = recenter, 
				left_atom_fixed = left_atom_fixed, flip_right = true, RNG = RNG)

			save(string("./states/", pot_type, "_state_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "state", state, "params", params)

			save(string("./data/", pot_type, "_period_avg_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "period_avg", period_avg)

			save(string("./data/", pot_type, "_fluxes_seed", seed, "_a", a, "_b", b, "_c", c, 
			"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
			"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
			"_period", recording_period, ".jld"), "p_flux", 
				p_flux, "e_flux", e_flux)

			save(string("./data/", pot_type, "_profiles_w_cov_seed", seed, "_a", a, "_b", b, "_c", c, 
			"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
			"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
			"_period", recording_period, ".jld"), "T_profile",
				T_profile, "cov_mat", cov_mat, "r_profile", r_profile, "e_profile", e_profile)
		else
			state, averages, T_profile, cov_mat, r_profile, e_profile, e_flux, p_flux = atom_chain(q_0, p_0, 
				grad_U, v, v_prime, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, T, N, burn_in_time, 
				save_cov_mat; rotar = rotar, pinning = alpha, recenter = recenter, 
				left_atom_fixed = left_atom_fixed, flip_right = true, RNG = RNG)

			save(string("./states/", pot_type, "_state_seed", seed, "_a", a, "_b", b, "_c", c, 
			"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
			"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
			"_period", recording_period, ".jld"), "state", state, "params", params)

			save(string("./data/", pot_type, "_fluxes_seed", seed, "_a", a, "_b", b, "_c", c, 
			"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
			"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
			"_period", recording_period, ".jld"), "p_flux", p_flux, "e_flux", e_flux)
	
			save(string("./data/", pot_type, "_profiles_w_cov_seed", seed, "_a", a, "_b", b, "_c", c, 
			"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
			"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
			"_period", recording_period, ".jld"), "T_profile", T_profile, 
				"cov_mat", cov_mat, "r_profile", r_profile, "e_profile", e_profile, "averages", 
				averages)
		end
			


	else
		if recording_period > 0
			state, period_avg, T_profile, q_var, r_profile, e_profile, e_flux, p_flux = atom_chain(q_0, p_0, 
				grad_U, v, v_prime, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, T, N, burn_in_time, 
				recording_period; rotar = rotar, pinning = alpha, recenter = recenter, 
				left_atom_fixed = left_atom_fixed, flip_right = true, RNG = RNG)

			save(string("./states/", pot_type, "_state_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "state", state, "params", params)

			save(string("./data/", pot_type, "_period_avg_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "period_avg", period_avg)

			save(string("./data/", pot_type, "_fluxes_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "p_flux", p_flux, "e_flux", e_flux)

			save(string("./data/", pot_type, "_profiles_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "T_profile", T_profile, "q_var", q_var, 
				"r_profile", r_profile, "e_profile", e_profile)

		else
			state, averages, T_profile, stretch_temp, tension_profile, q_var, r_profile, e_profile, e_flux, p_flux = atom_chain(q_0, p_0, 
				grad_U, v, v_prime, v_second, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, T, N, burn_in_time; 
				rotar = rotar, pinning = alpha, recenter = recenter, left_atom_fixed = left_atom_fixed, 
				flip_right = true, RNG = RNG)
			
			save(string("./states/", pot_type, "_state_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "state", state, "params", params)

			save(string("./data/", pot_type, "_fluxes_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "p_flux", p_flux, "e_flux", e_flux)

			save(string("./data/", pot_type, "_profiles_seed", seed, "_a", a, "_b", b, "_c", c, 
				"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
				"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
				"_period", recording_period, ".jld"), "T_profile", T_profile, "stretch_temp", 
				stretch_temp, "tension_profile", tension_profile, "q_var", q_var, "r_profile", 
				r_profile, "e_profile", e_profile, "averages", averages)

		end
	end

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
	tension = parse(Float64, ARGS[9])
	lambda = parse(Float64, ARGS[10])
	d = parse(Int64, ARGS[11])
	dt = parse(Float64, ARGS[12])
	T = parse(Float64, ARGS[13])
	burn_in_time = parse(Float64, ARGS[14])
	recording_period = parse(Float64, ARGS[15])
	save_cov_mat = parse(Bool, ARGS[16])
	seed = parse(Int, ARGS[17])


	print(string("Started: Potential = ", pot_type, ", a = ", a, ", b = ", b, ", c = ", c, 
		", alpha = ", alpha, ", omega = ", omega, ", tension = ", tension, ",\n",
		"T_l = ", T_l, ", gamma_l = ", gamma_l, ", lambda = ", lambda, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, " period = ", recording_period, 
		",\nSaving Covariance Matrix = ", save_cov_mat, ", seed = ", seed,
		"\n\n"))
	flush(stdout)


	main(pot_type, a, b, c, alpha, omega, d, T_l, gamma_l, tension, lambda, T, burn_in_time, dt, 
		recording_period, save_cov_mat, seed)
	print(string("Finished: Potential = ", pot_type, ", a = ", a, ", b = ", b, ", c = ", c, 
		", alpha = ", alpha, ", omega = ", omega, ", tension = ", tension, ",\n",
		"T_l = ", T_l, ", gamma_l = ", gamma_l, ", lambda = ", lambda, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, " period = ", recording_period,
		",\nSaving Covariance Matrix = ", save_cov_mat, ", seed = ", seed,
		"\n\n"))
	flush(stdout)
end


