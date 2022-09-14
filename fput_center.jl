using JLD

include("flux.jl")
include("atom_chains.jl")


function main(pot_type, a, b, c, alpha, omega, d, T_l, gamma_l, tension, lambda, T, burn_in_time, 
	time_step, skip, seed, cont, T_new = 0, burn_new = 0)
	gamma_r = 0.
	T_r = 0.

	F = get_forcing_cos(c, omega, d)

	left_atom_fixed = false
	recenter = false

	if pot_type == "FPUT"
		U, grad_U, v, v_prime, v_second = get_pot_fput(alpha, a, b)
		rotar = false

	elseif pot_type == "harmonic"
		U, grad_U, v, v_prime, v_second = get_pot_harm(alpha, a)
		rotar = false

	elseif pot_type == "cos"
		U, grad_U, v, v_prime, v_second = get_pot_cos(alpha, a)
		rotar = true

	else
		throw(ArgumentError(string("\"", pot_type, "\" not a valid potential type argument", 
            "\nValid potential type arguments: \"FPUT\", \"harmonic\", \"cos\"")))
	end


	if cont == 0
		Random.seed!(seed)
		RNG = copy(Random.default_rng())

		if pot_type == "FPUT"
			q_0 = reject_sampling_fput(RNG, d, v, a, b, T_l)
			p_0 = sqrt(T_l) * randn(d)
		elseif pot_type == "harmonic"
			q_0 = sqrt(T_l) * randn(d)
			p_0 = sqrt(T_l) * randn(d)
		elseif pot_type == "cos"
			q_0 = reject_sampling_cos(RNG, d, v, T_l)
			p_0 = sqrt(T_l) * randn(d)
		end

		T_new = T
		burn_new = burn_in_time
		old_params = ()
	else
		state = load(string("./states/", pot_type, "center_state_seed", seed, "_a", a, "_b", b, "_c", c, 
            "_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
            "_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
            "_skip", skip, "_cont", cont-1, ".jld"), "state")
		params = load(string("./states/", pot_type, "center_state_seed", seed, "_a", a, "_b", b, "_c", c, 
			"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
			"_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
			"_skip", skip, "_cont", cont-1, ".jld"), "params")

		old_params = load(string("./states/", pot_type, "center_state_seed", seed, "_a", a, "_b", b, 
			"_c", c, "_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", 
			lambda, "_tension", tension, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, 
            "_skip", skip, "_cont", cont-1, ".jld"), "old_params")

		old_params = (old_params..., params)
		
		q_0, p_0, RNG_state = state

		RNG = Xoshiro(RNG_state...)
	end

	params = Dict([("pot_type", pot_type), ("a", a), ("b", b), ("c", c), ("alpha", alpha), 
		("omega", omega), ("d", d), ("T_l", T_l), ("gamma_l", gamma_l), ("tension", tension), 
		("lambda", lambda), ("T", T_new), ("burn_in_time", burn_new), ("dt", dt), 
		("skip", skip), ("seed", seed), ("cont", cont)])

	N = floor(Int, T/time_step)


    state, center_mass, mv_avg_center_mass = atom_chain_center(q_0, p_0, grad_U, F, T_l, T_r, gamma_l, 
		gamma_r, tension, lambda, T_new, N, burn_new, skip; rotar, recenter, flip_right = true, 
		RNG = RNG)

    save(string("./data/", pot_type, "_center_seed", seed, "_a", a, "_b", b, "_c", c, "_alpha", alpha, 
		"_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_tension", tension, "_lambda", lambda, 
		"_d", d, "_T", T_new, "_burn", burn_new, "_dt", time_step, "_skip", skip, "_cont", cont, 
		".jld"), "center_mass", center_mass, "average_center", mv_avg_center_mass)

	save(string("./states/", pot_type, "center_state_seed", seed, "_a", a, "_b", b, "_c", c, 
		"_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, 
		"_tension", tension, "_d", d, "_T", T_new, "_burn", burn_new, "_dt", time_step, "_skip", 
		skip, "_cont", cont, ".jld"), "state", state, "params", params, "old_params", old_params)
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
	skip = parse(Int64, ARGS[15])
	seed = parse(Int, ARGS[16])
	cont = parse(Int, ARGS[17])
	T_new = parse(Int, ARGS[18])
	burn_new = parse(Int, ARGS[19])


	print(string("Started: Potential = ", pot_type, ", a = ", a, ", b = ", b, ", c = ", c, 
		", alpha = ", alpha, ", omega = ", omega, ", tension = ", tension, ",\n",
		"T_l = ", T_l, ", gamma_l = ", gamma_l, ", lambda = ", lambda, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, " skip = ", skip, 
		"\n\n"))
	flush(stdout)

    main(pot_type, a, b, c, alpha, omega, d, T_l, gamma_l, tension, lambda, T, burn_in_time, dt, 
		skip, seed, cont, T_new, burn_new)
    
	print(string("Finished: Potential = ", pot_type, ", a = ", a, ", b = ", b, ", c = ", c, 
		", alpha = ", alpha, ", omega = ", omega, ", tension = ", tension, ",\n",
		"T_l = ", T_l, ", gamma_l = ", gamma_l, ", lambda = ", lambda, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, " skip = ", skip,
		"\n\n"))
	flush(stdout)
end


