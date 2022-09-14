using JLD

include("flux.jl")
include("atom_chains.jl")


function main(pot_type, a, b, alpha, d, T_l, gamma_l, T_r, gamma_r, tension, lambda, T, 
	burn_in_time, time_step, skip, seed)

	F(t) = 0

	left_atom_fixed = false
	

	if pot_type == "FPUT"
		U, grad_U, v, v_prime, v_second = get_pot_fput(alpha, a, b)
		rotar = false
		recenter = (alpha == 0)

	elseif pot_type == "harmonic"
		U, grad_U, v, v_prime, v_second = get_pot_harm(alpha, a)
		rotar = false
		recenter = (alpha == 0)

	elseif pot_type == "cos"
		U, grad_U, v, v_prime, v_second = get_pot_cos(alpha, a)
		rotar = true

	else
		throw(ArgumentError(string("\"", pot_type, "\" not a valid potential type argument", 
            "\nValid potential type arguments: \"FPUT\", \"harmonic\", \"cos\"")))
	end


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


	N = floor(Int, T/time_step)


    trans_coef, running_trans_coef= atom_chain_trans_coef(q_0, p_0, grad_U, v_prime, F, T_l, T_r, 
    	gamma_l, gamma_r, tension, lambda, T, N, burn_in_time, skip; rotar, recenter, 
    	flip_right = false, RNG = RNG)

    save(string("./data/", pot_type, "_trans_coef_seed", seed, "_a", a, "_b", b, "_alpha", alpha, 
		"_Tl", T_l, "_Tr", T_r, "_gamma_l", gamma_l, "_gamma_r", gamma_r, "_tension", tension, "_lambda", 
		lambda, "_d", d, "_T", T, "_burn", burn_in_time, "_dt", time_step, "_skip", skip, ".jld"), 
    	"trans_coef", trans_coef, "running_trans_coef", running_trans_coef)
end


if abspath(PROGRAM_FILE) == @__FILE__
	pot_type = ARGS[1]
	a = parse(Float64, ARGS[2])
	b = parse(Float64, ARGS[3])
	alpha = parse(Float64, ARGS[4])
	T_l = parse(Float64, ARGS[5])
	T_r = parse(Float64, ARGS[6])
	gamma_l = parse(Float64, ARGS[7])
	gamma_r = parse(Float64, ARGS[8])
	tension = parse(Float64, ARGS[9])
	lambda = parse(Float64, ARGS[10])
	d = parse(Int64, ARGS[11])
	dt = parse(Float64, ARGS[12])
	T = parse(Float64, ARGS[13])
	burn_in_time = parse(Float64, ARGS[14])
	skip = parse(Int64, ARGS[15])
	seed = parse(Int, ARGS[16])


	print(string("Started: Potential = ", pot_type, ", a = ", a, ", b = ", b,", alpha = ", alpha, 
		", tension = ", tension, ",\n",
		"T_l = ", T_l, ", T_r = ", T_r, ", gamma_l = ", gamma_l, ", gamma_r =", gamma_r, 
		", lambda = ", lambda, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, " skip = ", skip, 
		"\n\n"))
	flush(stdout)

    main(pot_type, a, b, alpha, d, T_l, gamma_l, T_r, gamma_r, tension, lambda, T, 
		burn_in_time, dt, skip, seed)
    
	print(string("Finished: Potential = ", pot_type, ", a = ", a, ", b = ", b,", alpha = ", alpha, 
		", tension = ", tension, ",\n",
		"T_l = ", T_l, ", T_r = ", T_r, ", gamma_l = ", gamma_l, ", gamma_r =", gamma_r, 
		", lambda = ", lambda, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, " skip = ", skip, 
		"\n\n"))
	flush(stdout)
end


