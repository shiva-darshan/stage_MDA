using JLD

include("atom_chains.jl")


function main(pot_type, a, b, c, alpha, omega, d, T_l, gamma_l, lambda, T, burn_in_time, time_step, 
    period, q_upper, q_lower, q_dx, p_upper, p_lower, p_dx, r_upper, r_lower, r_dx)
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

	q_bins = q_lower:q_dx:q_upper
	p_bins = p_lower:p_dx:p_upper
	r_bins = r_lower:r_dx:r_upper

	N = floor(Int, T/time_step)

	if period == 0.
		q_histo, p_histo, r_histo = atom_chain_hist(q_0, p_0, grad_U, F, T_l, T_r, gamma_l, gamma_r, 
		lambda, T, N, burn_in_time, q_bins, p_bins, r_bins; rotar = rotar, recenter = recenter)

		save(string("./data/", pot_type, "_histos", "_qbins", q_bins, "_pbins", p_bins, "_rbins", r_bins, 
		"_a", a, "_b", b, "_c", c, "_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l,
		"_lambda", lambda, "_d", d, "_T", T, "_dt", time_step, ".jld"), "q_histo", q_histo, "p_histo",
		p_histo, "r_histo", r_histo)
	else
		q_histo, p_histo, r_histo = atom_chain_hist(q_0, p_0, grad_U, F, T_l, T_r, gamma_l, gamma_r, 
		lambda, T, N, burn_in_time, period, q_bins, p_bins, r_bins; rotar = rotar, recenter = recenter)

		save(string("./data/", pot_type, "_histos", "_qbins", q_bins, "_pbins", p_bins, "_rbins", r_bins, 
		"_a", a, "_b", b, "_c", c, "_alpha", alpha, "_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l,
		"_lambda", lambda, "_d", d, "_T", T, "_dt", time_step, "_period", period, ".jld"), "q_histo",
		q_histo, "p_histo", p_histo, "r_histo", r_histo)
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
	lambda = parse(Float64, ARGS[9])
	d = parse(Int64, ARGS[10])
	dt = parse(Float64, ARGS[11])
	T = parse(Float64, ARGS[12])
	burn_in_time = parse(Float64, ARGS[13])
	period = parse(Float64, ARGS[14])
	q_upper = parse(Float64, ARGS[15])
	q_lower = parse(Float64, ARGS[16])
	q_dx = parse(Float64, ARGS[17])
	p_upper = parse(Float64, ARGS[18])
	p_lower = parse(Float64, ARGS[19])
	p_dx = parse(Float64, ARGS[20])
	r_upper = parse(Float64, ARGS[21])
	r_lower = parse(Float64, ARGS[22])
	r_dx = parse(Float64, ARGS[23])
	


	print(string("Started: Trajs, Potential = ", pot_type, ", a = ", a, ", b = ", b, ", c = ", c, 
		", alpha = ", alpha, ", omega = ", omega, ",\n",
		"T_l = ", T_l, ", gamma_l = ", gamma_l, ", lambda = ", lambda, ",\n",
		"q_upper = ", q_upper, ", q_lower = ", q_lower, ", q_dx = ", q_dx, ",\n",
		"p_upper = ", p_upper, ", p_lower = ", p_lower, ", p_dx = ", p_dx, ",\n",
		"r_upper = ", r_upper, ", r_lower = ", r_lower, ", r_dx = ", r_dx, ",\n",
		"d = ", d, ", dt = ", dt, ", T = ", T, "\n\n"))
	flush(stdout)


	main(pot_type, a, b, c, alpha, omega, d, T_l, gamma_l, lambda, T, burn_in_time, dt, period, q_upper, 
		q_lower, q_dx, p_upper, p_lower, p_dx, r_upper, r_lower, r_dx)
	print(string("Finished: Trajs, Potential = ", pot_type, ", a = ", a, ", b = ", b, ", c = ", c, 
	", alpha = ", alpha, ", omega = ", omega, ",\n",
	"T_l = ", T_l, ", gamma_l = ", gamma_l, ", lambda = ", lambda, ",\n",
	"q_upper = ", q_upper, ", q_lower = ", q_lower, ", q_dx = ", q_dx, ",\n",
	"p_upper = ", p_upper, ", p_lower = ", p_lower, ", p_dx = ", p_dx, ",\n",
	"r_upper = ", r_upper, ", r_lower = ", r_lower, ", r_dx = ", r_dx, ",\n",
	"d = ", d, ", dt = ", dt, ", T = ", T, "\n\n"))
	flush(stdout)
end


