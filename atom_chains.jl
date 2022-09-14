using LinearAlgebra
using LinearAlgebra.BLAS
using Random
using Distributions
using LoopVectorization

function get_pot_fput(alpha::Float64, a::Float64, b::Float64, left_atom_fixed::Bool = false)
	"""
	Return FPUT with pinning potential and its gradient

	Parameters
	-----------

	alpha : intensity of pinning potential

	a : cubic coefficient of FPUT potential

	b : quartic coefficient of FPUT potential

	"""


	#FPUT potential local interactions
	v(r) = r^2/2 + a*r^3/3 + b*r^4/4
	v_prime(r) = r + a*r^2 + b*r^3 
	v_second(r) = 1 + 2*a*r + 3*b*r^2

	function U_fput(x)
		pinning = alpha * (x'*x)/2
		interaction = left_atom_fixed * v(x[1]) + sum(v.(x[2:end] - x[1:end-1]))
		return pinning + interaction
	end


	function grad_U_fput!(pot::Array{Float64, 1}, left_interaction::Array{Float64, 1}, 
		right_interaction::Array{Float64, 1}, x::Array{Float64, 1})
		pot .= alpha .* x # gradient of confining potential
		
		#gradient of interaction potential
		left_interaction[1] = left_atom_fixed * v_prime(x[1])
		right_interaction[1] = -v_prime(x[2] - x[1])
		@avx for i = 2:(length(x)-1)
			left_interaction[i] = v_prime(x[i] - x[i-1])
			right_interaction[i] = -v_prime(x[i+1] - x[i])
		end
		left_interaction[end] = v_prime(x[end]-x[end-1])
		pot .+= left_interaction
		pot .+= right_interaction
		return pot
	end

	return U_fput, grad_U_fput!, v, v_prime, v_second

end

function reject_sampling_fput(RNG, d, v, a, b, T_mid)
	#only works for a = -2.5, b = 1
	#TODO update to work for any parameter choice

	M = 25
	M = max(M/T_mid, M * T_mid)
	xs = LinRange(-10, 10, 20001)
	Z = sum(exp.(-v.(xs)/T_mid))*(xs[2] - xs[1])
	m = (-a + sign(a/(-2.5)) * sqrt(a^2 - 4b))/(4*b)
	norm_pdf(x) = exp(-(x -m)^2/2)/sqrt(2 * pi)

	#density of FPUT bond diameter probability
	f(x) = exp(-v(x)/T_mid)/Z

	out = zeros(d)
	samples = 0
	while samples < d
		prop = 1 + randn(RNG)
		U = rand(RNG)

		if U < f(prop)/(M * norm_pdf(prop))
			samples += 1
			out[samples] = prop
		end
	end

	out = cumsum(out)

	return out
end

function get_pot_harm(alpha, a, left_atom_fixed = false)
	"""
	Return harmonic with pinning potential and its gradient

	Parameters
	-----------

	alpha : intensity of pinning potential

	a : coefficient of harmonic potential

	"""

	v(r) = a * r^2/2
	v_prime(r) = a * r
	v_second(r) = a

	function U_harm(x)
		pinning = alpha * (x'*x)/2
		interaction = (a/2) * (left_atom_fixed * dot(x[1], x[1]) + dot(x[2:end] - x[1:end-1],x[2:end] - x[1:end-1]))
		return pinning + interaction
	end
 

	function grad_U_harm!(pot::Array{Float64, 1}, left_interaction::Array{Float64, 1}, 
		right_interaction::Array{Float64, 1}, x::Array{Float64, 1})
		pot .= alpha .* x # gradient of confining potential
		
		#gradient of interaction potential
		left_interaction[1] = a * left_atom_fixed * x[1]
		right_interaction[1] = a * (x[1] - x[2]) #swapped x[1] and x[2] since x->x is odd
		@avx for i = 2:(length(x)-1)
			left_interaction[i] = a * (x[i] - x[i-1])
			right_interaction[i] = a * (x[i] - x[i+1]) #swapped x[i] and x[i+1] since x->x is odd
		end
		left_interaction[end] = a * (x[end]-x[end-1])
		pot .+= left_interaction
		pot .+= right_interaction
		return pot
	end

	return U_harm, grad_U_harm!, v, v_prime, v_second

end

function get_pot_cos(alpha::Float64, a::Float64, left_atom_fixed::Bool = false)
	"""
	Return cos potential and its gradient

	not pinned on the left hand side

	Parameters
	-----------

	"""

	v(r) = a * (1 - cos(r))
	v_prime(r) = a * sin(r)
	v_second(r) = a * cos(r)

	function U_cos(x)
		interaction = a*sum(v.(x[2:end] - x[1:end-1])) + left_atom_fixed*a*v(x[1])
		pinning =  alpha*x'*x/2
		return interaction + pinning
	end


	function grad_U_cos!(pot::Array{Float64, 1}, left_interaction::Array{Float64, 1}, 
		right_interaction::Array{Float64, 1}, x::Array{Float64, 1})
		pot .= alpha .* x # gradient of confining potential
		
		#gradient of interaction potential
		left_interaction[1] = left_atom_fixed * v_prime(x[1])
		right_interaction[1] = v_prime(x[1] - x[2]) #swapped x[1] and x[2] since sin is odd
		@avx for i = 2:(length(x)-1)
			left_interaction[i] = v_prime(x[i] - x[i-1])
			right_interaction[i] = v_prime(x[i] - x[i+1]) #swapped x[i] and x[i+1] since sin is odd
		end
		left_interaction[end] = v_prime(x[end]-x[end-1])
		pot .+= left_interaction
		pot .+= right_interaction
		return pot
	end

	return U_cos, grad_U_cos!, v, v_prime, v_second

end

function get_forcing_cos(c::Float64, omega::Float64, d::Int64)
	function F_cos(t::Float64)
		return (c/sqrt(d)) * cos(2 * pi * omega * t)
	end

	return F_cos
end

function reject_sampling_cos(RNG, d, v, T_mid)
	#only works for a = -2.5, b = 1
	#TODO update to work for any parameter choice

	M = 2.5
	M = max(M/T_mid, M * T_mid)
	xs = LinRange(0, 2 * pi, 10001)
	Z = sum(exp.(-v.(xs)/T_mid))*(xs[2] - xs[1])


	#density of cosine bond diameter probability
	f(x) = exp(-v(x)/T_mid)/Z


	out = zeros(d)
	samples = 0
	while samples < d
		prop = 2 * pi * rand(RNG)
		U = rand(RNG)

		if U < f(prop)/(M/(2 * pi))
			samples += 1
			out[samples] = prop
		end
	end

	out = cumsum(out)

	out = ((out .% (2*pi)) .+ (2*pi)) .% (2*pi)
	return out
end


function integration_step!(q_n::Array{Float64,1}, p_n::Array{Float64,1}, 
	ham_part::Array{Float64,1}, grad_pot::Array{Float64,1}, 
	left_interaction::Array{Float64,1}, right_interaction::Array{Float64,1}, i::Int64, 
	grad_U!::Function, F::Function, tension::AbstractFloat, T_l::AbstractFloat, T_r::AbstractFloat, 
	gamma_l::AbstractFloat, gamma_r::AbstractFloat, dt::AbstractFloat; RNG::AbstractRNG)
	"""

	One step of Velocity Verlet integration of atom chain 
	"""
	q_n .+= dt/2 .* p_n

	ham_part .= .-grad_U!(grad_pot, left_interaction, right_interaction, q_n) .* dt


	p_n[2:end-1] .+= @view(ham_part[2:end-1])


	#thermostated boundaries
	p_n[1] += ham_part[1] - gamma_l*p_n[1]*dt + sqrt(2*T_l*gamma_l*dt)*randn(RNG) 
	p_n[end] += ham_part[end] - gamma_r*p_n[end]*dt + sqrt(2*T_r*gamma_r*dt)*randn(RNG) 

	#forcing at right endpoint
	p_n[end] += F(dt * i) * dt + tension * dt


	q_n .+= dt/2 .* p_n

	return nothing
end


function update_fluxes!(energy_flux, momentum_flux, i, r, p_n, v_prime::Function, F::Function, 
	gamma_l::AbstractFloat, gamma_r::AbstractFloat, T_l::AbstractFloat, T_r::AbstractFloat, 
	d::Int64, dt::AbstractFloat)
	"""

	"""
	energy_flux[1] += gamma_l *(T_l - p_n[1]^2)
	momentum_flux[1] += -gamma_l*p_n[1]

	energy_flux[end] += -(1/sqrt(d)) * F(dt * i) * p_n[end] - gamma_r*(T_r - p_n[end]^2)
	momentum_flux[end] += gamma_r*p_n[end] - (1/sqrt(d)) * F(dt * i) 

	@avx for j = 2:d
		energy_flux[j] += -p_n[j-1]*v_prime(r[j-1])
		momentum_flux[j] += -v_prime(r[j-1])
	end

	return nothing
end


function update_profiles!(energy_profile, momentum_profile, stretch_profile, q_profile, q2_stat, p2_stat,
	q_n, p_n, r_n, v::Function, pinning::AbstractFloat, d::Int64, left_atom_fixed::Bool, recenter::Bool)
	"""

	"""
	energy_profile[1] += !recenter*left_atom_fixed*v(q_n[1]) + !recenter*pinning*q_n[1]^2/2 + p_n[1]^2/2
	@avx for j = 2:d
		energy_profile[j] += v(r_n[j-1]) + pinning*q_n[j]^2/2 + p_n[j]^2/2 
	end

	stretch_profile .+= r_n
	
	q_profile .+= q_n
	momentum_profile .+= p_n
	q2_stat .+= q_n.^2
	p2_stat .+= p_n.^2

	return nothing
end


function update_profiles!(energy_profile, momentum_profile, stretch_profile, q_profile, q2_stat, 
	p2_stat, v_first, v_first_sq, v_second, q_n, p_n, r_n, v::Function, v_p::Function, v_s::Function, 
	pinning::AbstractFloat, d::Int64, left_atom_fixed::Bool, recenter::Bool)
	"""

	"""
	energy_profile[1] += !recenter*left_atom_fixed*v(q_n[1]) + !recenter*pinning*q_n[1]^2/2 + p_n[1]^2/2
	@avx for j = 2:d
		energy_profile[j] += v(r_n[j-1]) + pinning*q_n[j]^2/2 + p_n[j]^2/2
		v_first[j-1] += v_p(r_n[j-1])
		v_first_sq[j-1] += v_p(r_n[j-1])^2
		v_second[j-1] += v_s(r_n[j-1])
	end

	stretch_profile .+= r_n
	
	q_profile .+= q_n
	momentum_profile .+= p_n
	q2_stat .+= q_n.^2
	p2_stat .+= p_n.^2

	return nothing
end

function update_profiles!(energy_profile, momentum_profile, stretch_profile, q_profile, p2_stat, q_n,
	p_n, r_n, v::Function, pinning::AbstractFloat, d::Int64, left_atom_fixed::Bool, recenter::Bool)
	"""

	"""
	energy_profile[1] += !recenter*left_atom_fixed*v(q_n[1]) + !recenter*pinning*q_n[1]^2/2 + p_n[1]^2/2
	@avx for j = 2:d
		energy_profile[j] += v(r_n[j-1]) + pinning*q_n[j]^2/2 + p_n[j]^2/2 
	end

	stretch_profile .+= r_n
	
	q_profile .+= q_n
	momentum_profile .+= p_n
	p2_stat .+= p_n.^2

	return nothing
end

function update_profiles!(momentum_profile, p2_stat, bulk_energy_flux, i, r, p_n, v_prime::Function, d::Int64)
	"""

	"""

	@avx for j = 1:(d-1)
		bulk_energy_flux[j] += -p_n[j]*v_prime(r[j])
	end

	momentum_profile .+= p_n
	p2_stat .+= p_n.^2

	return nothing
end

@inline function update_counts!(counts, bins, dx, x, right_open::Bool = true)
	"""
	"""
	ind = max(floor(Int64, (x - bins[1])/dx) + right_open, 0) + 1

	if ind < length(counts) 
		counts[ind] += 1
	else
		counts[end] += 1
	end
	
	return nothing
end

function atom_chain(q_0::Array, p_0::Array, grad_U!::Function, v::Function, v_prime::Function,
	F::Function, T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat,
	tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, N::Integer, burn_in_time::AbstractFloat; 
	rotar::Bool = false, pinning::AbstractFloat = 0., recenter::Bool = true, left_atom_fixed::Bool = !recenter, 
	flip_left::Bool = false, flip_right::Bool = false, RNG::AbstractRNG = copy(Random.default_rng()))

	"""

	"""


	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	r_n = zeros(d-1) #bond distance
	
	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)

	temp_profile = zeros(d)
	energy_profile = zeros(d)
	q_profile = zeros(d)
	momentum_profile = zeros(d)
	stretch_profile = zeros(d-1)
	p2_stat = zeros(d)
	q2_stat = zeros(d)

	energy_flux = zeros(d+1)
	momentum_flux = zeros(d+1)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Array{Int64, 1}()
	#burn in 
	for i in 1:K
		#momentum flips
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
		grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)

	end

	for i in 1:N

		#momentum flips
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)


		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])


		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		update_profiles!(energy_profile, momentum_profile, stretch_profile, q_profile, q2_stat, 
			p2_stat, q_n, p_n, r_n, v, pinning, d, left_atom_fixed, recenter)
		update_fluxes!(energy_flux, momentum_flux, i, r_n, p_n, v_prime, F, gamma_l, gamma_r, T_l, 
			T_r, d, dt)
	end

	RNG_state = ntuple(fieldcount(typeof(RNG))) do i
		getfield(rng, i)
	end

	state = (q_n, p_n, RNG_state)

	energy_profile /= N
	momentum_profile /= N
	stretch_profile /= N
	q_profile /= N
	temp_profile = p2_stat/N - momentum_profile.^2
	q_var = q2_stat/N - q_profile.^2

	energy_flux /= N
	momentum_flux /= N

	averages = zeros(d, 2)
	averages[:, 1] = q_profile
	averages[:, 2] = momentum_profile

	return state, averages, temp_profile, q_var, stretch_profile, energy_profile, energy_flux, momentum_flux
end


function atom_chain(q_0::Array, p_0::Array, grad_U!::Function, v::Function, v_prime::Function, 
	v_second::Function, F::Function, T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, 
	gamma_r::AbstractFloat, tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, 
	N::Integer, burn_in_time::AbstractFloat; rotar::Bool = false, pinning::AbstractFloat = 0., 
	recenter::Bool = true, left_atom_fixed::Bool = !recenter, flip_left::Bool = false, 
	flip_right::Bool = false, RNG::AbstractRNG = copy(Random.default_rng()))

	"""

	"""


	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	r_n = zeros(d-1) #bond distance
	
	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)

	temp_profile = zeros(d)
	energy_profile = zeros(d)
	q_profile = zeros(d)
	momentum_profile = zeros(d)
	stretch_profile = zeros(d-1)
	v_p = zeros(d-1)
	v_p_sq = zeros(d-1)
	v_pp = zeros(d-1)
	p2_stat = zeros(d)
	q2_stat = zeros(d)

	energy_flux = zeros(d+1)
	momentum_flux = zeros(d+1)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Array{Int64, 1}()
	#burn in 
	for i in 1:K
		#momentum flips
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
		grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)

	end

	for i in 1:N

		#momentum flips
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)


		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])


		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		update_profiles!(energy_profile, momentum_profile, stretch_profile, q_profile, q2_stat, 
			p2_stat, v_p, v_p_sq, v_pp, q_n, p_n, r_n, v, v_prime, v_second, pinning, d, 
			left_atom_fixed, recenter)

		update_fluxes!(energy_flux, momentum_flux, i, r_n, p_n, v_prime, F, gamma_l, gamma_r, T_l, 
			T_r, d, dt)
	end

	RNG_state = ntuple(fieldcount(typeof(RNG))) do i
		getfield(RNG, i)
	end

	state = (q_n, p_n, RNG_state)

	energy_profile /= N
	momentum_profile /= N
	stretch_profile /= N
	q_profile /= N
	temp_profile = p2_stat/N - momentum_profile.^2
	stretch_temp = (v_p_sq .- (v_p.^2)./N) ./ v_pp
	tau = v_p ./ N
	q_var = q2_stat/N - q_profile.^2

	energy_flux /= N
	momentum_flux /= N

	averages = zeros(d, 2)
	averages[:, 1] = q_profile
	averages[:, 2] = momentum_profile

	return state, averages, temp_profile, stretch_temp, tau, q_var, stretch_profile, energy_profile, energy_flux, momentum_flux
end

function atom_chain(q_0::Array, p_0::Array, grad_U!::Function, v::Function, v_prime::Function,
	F::Function, T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat,
	tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, N::Integer, burn_in_time::AbstractFloat, 
	save_S_mat::Bool; rotar::Bool = false, pinning::AbstractFloat = 0., recenter::Bool = true, 
	left_atom_fixed::Bool = !recenter, flip_left::Bool = false, flip_right::Bool = false,
	RNG::AbstractRNG = copy(Random.default_rng()))

	"""

	"""


	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	x_n = zeros(2d) #state vector
	r_n = zeros(d-1) #bond lengths

	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)

	temp_profile = zeros(d)
	energy_profile = zeros(d)
	q_profile = zeros(d)
	momentum_profile = zeros(d)
	stretch_profile = zeros(d-1)
	p2_stat = zeros(d)

	cross_moment = zeros(2d, 2d)

	energy_flux = zeros(d+1)
	momentum_flux = zeros(d+1)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end


	signs_to_change = Vector{Int64}()

	#burn in 
	for i in 1:K
		#momentum flip
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
		grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)
	end

	for i in 1:N
		#momentum flip
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)

		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])

		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		x_n[1:d] = q_n
		x_n[d+1:end] = p_n
		syr!('U', 1., x_n, cross_moment)


		update_profiles!(energy_profile, momentum_profile, stretch_profile, q_profile, p2_stat, q_n, 
		p_n, r_n, v, pinning, d, left_atom_fixed, recenter)
		update_fluxes!(energy_flux, momentum_flux, i, r_n, p_n, v_prime, F, gamma_l, gamma_r, T_l, T_r, d, dt)
	end

	RNG_state = ntuple(fieldcount(typeof(RNG))) do i
		getfield(RNG, i)
	end

	state = (q_n, p_n, RNG_state)

	
	stretch_profile /= N
	energy_profile /= N
	q_profile /= N
	momentum_profile /= N
	temp_profile = p2_stat/N - momentum_profile.^2
	cross_moment = Symmetric(cross_moment ./ N) 
	S_mat = cross_moment - Symmetric([q_profile; momentum_profile] * [q_profile; momentum_profile]')
	energy_flux /= N
	momentum_flux /= N

	averages = zeros(d, 2)
	averages[:, 1] = q_profile
	averages[:, 2] = momentum_profile

	return state, averages, temp_profile, S_mat, stretch_profile, energy_profile, energy_flux, momentum_flux
end


function atom_chain(q_0::Array, p_0::Array, grad_U!::Function, v::Function, v_prime::Function,
	F::Function, T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat,
	tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, N::Integer, burn_in_time::AbstractFloat,
	period::AbstractFloat; rotar::Bool = false, pinning::AbstractFloat = 0., recenter::Bool = true, 
	left_atom_fixed::Bool = !recenter, flip_left::Bool = false, flip_right::Bool = false, 
	RNG::AbstractRNG = copy(Random.default_rng()))

	"""


	"""

	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	r_n = zeros(d-1) #bond lengths

	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)
	M = floor(Int, period/dt)

	temp_profile = zeros(d, M)
	energy_profile = zeros(d, M)
	q_avg = zeros(d, M)
	p_avg = zeros(d, M)
	stretch_profile = zeros(d-1, M)
	p2_stat = zeros(d, M)
	q2_stat = zeros(d, M)

	energy_flux = zeros(d+1, M)
	momentum_flux = zeros(d+1, M)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Vector{Int64}()

	#burn in 
	for i in 1:K
		#momentum flip
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, T_l, T_r, tension, gamma_l, gamma_r, dt; RNG = RNG)
	end

	for i in 1:N
		#momentum flip
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, T_l, T_r, tension, gamma_l, gamma_r, dt; RNG = RNG)


		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])

		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end


		update_profiles!(view(energy_profile, :, (i-1)%M + 1), view(p_avg, :, (i-1)%M + 1), 
		view(stretch_profile, :, (i-1)%M + 1), view(q_avg, :, (i-1)%M + 1), view(q2_stat, :, (i-1)%M + 1),
		view(p2_stat, :, (i-1)%M + 1), q_n, p_n, r_n, v, pinning, d, left_atom_fixed, recenter)
		update_fluxes!(view(energy_flux, :, (i-1)%M + 1), view(momentum_flux, :, (i-1)%M + 1), i, 
		r_n, p_n, v_prime, F, gamma_l, gamma_r, T_l, T_r, d, dt)
	end

	RNG_state = ntuple(fieldcount(typeof(RNG))) do i
		getfield(RNG, i)
	end

	state = (q_n, p_n, RNG_state)
	
	stretch_profile /= (N/M)
	energy_profile /= (N/M) 
	p_avg /= (N/M)
	q_avg /= (N/M)
	temp_profile = p2_stat/(N/M) - p_avg.^2
	q_var = q2_stat/(N/M) - q_avg.^2

	period_avg = zeros(d, M, 2)
	period_avg[:, :, 1] = q_avg
	period_avg[:, :, 2] = p_avg

	energy_flux /= (N/M)
	momentum_flux /= (N/M)

	return state, period_avg, temp_profile, q_var, stretch_profile, energy_profile, energy_flux, momentum_flux
end


function atom_chain(q_0::Array, p_0::Array, grad_U!::Function, v::Function, v_prime::Function,
	F::Function, T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat,
	tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, N::Integer, burn_in_time::AbstractFloat,
	period::AbstractFloat, save_S_mat::Bool; rotar::Bool = false, pinning::AbstractFloat = 0., 
	recenter::Bool = true, left_atom_fixed::Bool = !recenter, flip_left::Bool = false, 
	flip_right::Bool = false, RNG::AbstractRNG = copy(Random.default_rng()))

	"""


	"""

	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	x_n = zeros(2d) #state vector
	r_n = zeros(d-1) #bond lengths

	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)
	M = floor(Int, period/dt)

	temp_profile = zeros(d, M)
	energy_profile = zeros(d, M)
	q_avg = zeros(d, M)
	p_avg = zeros(d, M)
	stretch_profile = zeros(d-1, M)
	p2_stat = zeros(d, M)
	S_mat = zeros(2d, 2d, M)

	energy_flux = zeros(d+1, M)
	momentum_flux = zeros(d+1, M)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Vector{Int64}()

	#burn in 
	for i in 1:K
		#momentum flip
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG)
	end

	for i in 1:N
		#momentum flip
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG)


		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])

		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		x_n[1:d] = q_n
		x_n[d+1:end] = p_n
		syr!('U', 1., x_n, @view S_mat[ :, :, (i-1)%M + 1])


		update_profiles!(view(energy_profile, :, (i-1)%M + 1), view(p_avg, :, (i-1)%M + 1), 
		view(stretch_profile, :, (i-1)%M + 1), view(q_avg, :, (i-1)%M + 1),
		view(p2_stat, :, (i-1)%M + 1), q_n, p_n, r_n, v, pinning, d, left_atom_fixed, recenter)
		update_fluxes!(view(energy_flux, :, (i-1)%M + 1), view(momentum_flux, :, (i-1)%M + 1), i, 
		r_n, p_n, v_prime, F, gamma_l, gamma_r, T_l, T_r, d, dt)
	end

	RNG_state = ntuple(fieldcount(typeof(RNG))) do i
		getfield(RNG, i)
	end

	state = (q_n, p_n, RNG_state)
	
	stretch_profile /= (N/M)
	energy_profile /= (N/M) 
	p_avg /= (N/M)
	q_avg /= (N/M)
	temp_profile = p2_stat/(N/M) - p_avg.^2
	x_profile = zeros(2d, M)
	x_profile[1:d, :] = q_avg
	x_profile[d+1:end, :] = p_avg

	energy_flux /= (N/M)
	momentum_flux /= (N/M)

	S_mat ./= (N/M)

	for j = 1:M 
		syr!('U', -1., view(x_profile, :, j), view(S_mat, :, :, j))
	end

	period_avg = zeros(d, M, 2)
	period_avg[:, :, 1] = q_avg
	period_avg[:, :, 2] = p_avg

	return state, period_avg, temp_profile, S_mat, stretch_profile, energy_profile, energy_flux, momentum_flux
end


function atom_chain_center(q_0::Array, p_0::Array, grad_U!::Function, F::Function, T_l::AbstractFloat, 
	T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat, tension::AbstractFloat, 
	lambda::AbstractFloat, T::AbstractFloat, N::Integer, burn_in_time::AbstractFloat, 
	skip::Integer; rotar::Bool = false, recenter::Bool = true, flip_left::Bool = false, 
	flip_right::Bool = false, RNG::AbstractRNG = copy(Random.default_rng()))

	"""

	"""


	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	r_n = zeros(d-1) #bond distance
	
	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)

	baricentre_q = zeros(N ÷ skip)
	running_avg_q = zeros(N ÷ skip)
	cumsum_q = 0.

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Array{Int64, 1}()
	#burn in 
	for i in 1:K
		#momentum flips
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
		grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)

		q_n .-= q_n[1] #shift 1st atom back to zero

	end

	for i in 1:N

		#momentum flips
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)


		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])


		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		cumsum_q += mean(q_n)

		if i % skip == 0
			baricentre_q[i÷skip] = mean(q_n)
			running_avg_q[i÷skip] = cumsum_q/i
		end
	end

	RNG_state = ntuple(fieldcount(typeof(RNG))) do i
		getfield(RNG, i)
	end

	state = (q_n, p_n, RNG_state)

	return state, baricentre_q, running_avg_q
end


function atom_chain_trans_coef(q_0::Array, p_0::Array, grad_U!::Function, v_prime::Function, F::Function, 
	T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat, 
	tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, N::Integer, 
	burn_in_time::AbstractFloat, skip::Integer; rotar::Bool = false, recenter::Bool = true, 
	flip_left::Bool = false, flip_right::Bool = false, RNG::AbstractRNG = copy(Random.default_rng()))

	"""

	"""


	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	r_n = zeros(d-1) #bond distance
	
	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)

	trans_coef = zeros(d-1)
	running_trans_coef = zeros(d-1, N ÷ skip)
	bulk_e_flux = zeros(d-1)
	p2_stat = zeros(d)
	momentum_profile = zeros(d)
	temp_profile = zeros(d)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Array{Int64, 1}()
	#burn in 
	for i in 1:K
		#momentum flips
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
		grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)

		q_n .-= q_n[1] #shift 1st atom back to zero

	end

	for i in 1:N

		#momentum flips
		randsubseq!(RNG, signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt; RNG = RNG)


		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])


		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		update_profiles!(momentum_profile, p2_stat, bulk_e_flux, i, r_n, p_n, v_prime::Function, d::Int64)


		# trans_coef .= mean(bulk_e_flux) ./ ((@view p2_stat[2:end]) .- (@view p2_stat[1:end-1]) .- 
		# 	(@view momentum_profile[2:end]).^2 .+ (@view momentum_profile[1:end-1]).^2)


		if i % skip == 0
			temp_profile .= p2_stat ./i .- momentum_profile.^2 ./i
			J = mean(bulk_e_flux)/i

			running_trans_coef[:, i÷skip] .= J ./ ((@view temp_profile[2:end]) .-
			 (@view temp_profile[1:end-1]))
		end
	end

	temp_profile .= p2_stat./N .- momentum_profile.^2 ./N
	J = mean(bulk_e_flux)/N

	trans_coef .= J ./ ((@view temp_profile[2:end]) .- (@view temp_profile[1:end-1]))


	return trans_coef, running_trans_coef
end


function atom_chain_traj(q_0::Array, p_0::Array, grad_U!::Function,
	F::Function, T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat,
	tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, N::Integer, burn_in_time::AbstractFloat; 
	rotar::Bool = false, recenter::Bool = true, flip_left::Bool = false, flip_right::Bool = true)

	"""

	"""


	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	
	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)

	q_traj = zeros(d, N)
	p_traj = zeros(d, N)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Array{Int64, 1}()
	#burn in 
	for i in 1:K
		#momentum flips
		randsubseq!(signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt)
	end

	for i in 1:N
		#momentum flips
		randsubseq!(signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt)

		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		q_traj[:, i] .= q_n
		p_traj[:, i] .= p_n
	end

	return q_traj, p_traj
end


function atom_chain_hist(q_0::Array, p_0::Array, grad_U!::Function,
	F::Function, T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat,
	tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, N::Integer, burn_in_time::AbstractFloat, 
	q_bins::StepRangeLen, p_bins::StepRangeLen, r_bins::StepRangeLen; rotar::Bool = false, 
	recenter::Bool = true, flip_left::Bool = false, flip_right::Bool = true, return_counts::Bool = false)

	"""

	"""
	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	r_n = zeros(d-1) #bond lengths

	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)

	q_counts = zeros(Int64, length(q_bins) + 1, d)
	q_dx = Float64(q_bins.step)
	p_counts = zeros(Int64, length(p_bins) + 1, d)
	p_dx = Float64(p_bins.step)
	r_counts = zeros(Int64, length(r_bins) + 1, d-1)
	r_dx = Float64(r_bins.step)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Vector{Int64}()

	#burn in 
	for i in 1:K
		#momentum flip
		randsubseq!(signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt)
	end

	for i in 1:N
		#momentum flip
		randsubseq!(signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt)

		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])

		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		for j = 1:d-1
			update_counts!(view(q_counts, :, j), q_bins, q_dx, q_n[j])
			update_counts!(view(p_counts, :, j), p_bins, p_dx, p_n[j])
			update_counts!(view(r_counts, :, j), r_bins, r_dx, r_n[j])
		end

		update_counts!(view(q_counts, :, d), q_bins, q_dx, q_n[end])
		update_counts!(view(p_counts, :, d), p_bins, p_dx, p_n[end])
	end

	q_histo = q_counts/(N*q_dx)
	p_histo = p_counts/(N*p_dx)
	r_histo = r_counts/(N*r_dx)

	if return_counts
		return q_histo, p_histo, r_histo, q_counts, p_counts, r_counts
	else
		return q_histo, p_histo, r_histo
	end
end

function atom_chain_hist(q_0::Array, p_0::Array, grad_U!::Function,
	F::Function, T_l::AbstractFloat, T_r::AbstractFloat, gamma_l::AbstractFloat, gamma_r::AbstractFloat,
	tension::AbstractFloat, lambda::AbstractFloat, T::AbstractFloat, N::Integer, burn_in_time::AbstractFloat, 
	period::AbstractFloat, q_bins::StepRangeLen, p_bins::StepRangeLen, r_bins::StepRangeLen; rotar::Bool = false, 
	recenter::Bool = true, flip_left::Bool = false, flip_right::Bool = true, return_counts::Bool = false)

	"""

	"""
	d = length(q_0)
	q_n = copy(q_0) #position
	p_n = copy(p_0) #momentum
	r_n = zeros(d-1) #bond lengths

	ham_part = zeros(d)
	grad_pot = zeros(d)
	left_interaction = zeros(d)
	right_interaction = zeros(d)

	dt = T/N
	N = floor(Int, T/dt)
	dt = T/N
	K = floor(Int, burn_in_time/dt)
	M = floor(Int, period/dt)

	q_counts = zeros(Int64, M, d, length(q_bins) + 1)
	q_dx = Float64(q_bins.step)
	p_counts = zeros(Int64, M, d, length(p_bins) + 1)
	p_dx = Float64(p_bins.step)
	r_counts = zeros(Int64, M, d-1, length(r_bins) + 1)
	r_dx = Float64(r_bins.step)

	p_sign_change = - 1/2 * expm1(-2*lambda*dt)

	#fix which atoms will be flipped
	flippable_inds = 2:(d-1)
	if flip_left
		flippable_inds = vcat(1, flippable_inds)
	end
	if flip_right
	    flippable_inds = vcat(flippable_inds, d)
	end

	signs_to_change = Vector{Int64}()

	#burn in 
	for i in 1:K
		#momentum flip
		randsubseq!(signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] *= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt)
	end

	for i in 1:N
		#momentum flip
		randsubseq!(signs_to_change, flippable_inds, p_sign_change) #which atoms had a momentum sign change
		p_n[signs_to_change] .*= -1

		integration_step!(q_n, p_n, ham_part, grad_pot, left_interaction, right_interaction, i, 
			grad_U!, F, tension, T_l, T_r, gamma_l, gamma_r, dt)

		#update bond distance
		r_n .= (@view q_n[2:end]) .- (@view q_n[1:end-1])

		if rotar 
			#send back to torus if rotar chain
			q_n .= ((q_n .% (2*pi)) .+ (2*pi)) .% (2*pi)
		end

		if recenter
			q_n .-= q_n[1] #shift 1st atom back to zero
		end

		for j = 1:d-1
			update_counts!(view(q_counts, (i -1)%M + 1, j, :), q_bins, q_dx, q_n[j])
			update_counts!(view(p_counts, (i -1)%M + 1, j, :), p_bins, p_dx, p_n[j])
			update_counts!(view(r_counts, (i -1)%M + 1, j, :), r_bins, r_dx, r_n[j])
		end

		update_counts!(view(q_counts, (i -1)%M + 1, d, :), q_bins, q_dx, q_n[end])
		update_counts!(view(p_counts, (i -1)%M + 1, d, :), p_bins, p_dx, p_n[end])
	end

	q_histo = q_counts/((N/M)*q_dx)
	p_histo = p_counts/((N/M)*p_dx)
	r_histo = r_counts/((N/M)*r_dx)

	if return_counts
		return q_histo, p_histo, r_histo, q_counts, p_counts, r_counts
	else
		return q_histo, p_histo, r_histo
	end
end