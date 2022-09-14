using LinearAlgebra
using Random
using Distributions


function compute_momentum_flux(traj_q, traj_p, v_prime, offset)
	r = hcat(traj_q[offset:end, 1] ,traj_q[offset:end, 2:end] - traj_q[offset:end, 1:end-1])

	tau = size(r)[1]

	summand = - v_prime.(r)

	return (sum(summand, dims = 1)/tau)[1, :]

end

function compute_energy_flux(traj_q, traj_p, v_prime, offset)

	#bond lengths
	r = hcat(traj_q[offset:end, 1] ,traj_q[offset:end, 2:end] - traj_q[offset:end, 1:end-1])

	p = traj_p[offset:end, :]

	tau = size(p)[1]

	summand = - p .* v_prime.(r)

	return (sum(summand, dims = 1)/tau)[1, :]
end

function compute_temp_profile(traj_p, offset)

	p = traj_p[offset:end, :]

	tau = size(p)[1]

	p_stat = sum(p, dims = 1)/tau

	p_sq_stat = sum(p.^2, dims = 1)/tau

	return (p_sq_stat - (p_stat).^2)[1, :]

end

function compute_momentum_profile(traj_p, offset)

	p = traj_p[offset:end, :]

	tau = size(p)[1]

	p_stat = sum(p, dims = 1)/tau

	return p_stat[1, :]

end

function compute_energy_profile(traj_q, traj_p, v, offset)

	r = hcat(traj_q[offset:end, 1] ,traj_q[offset:end, 2:end] - traj_q[offset:end, 1:end-1])

	p = traj_p[offset:end, :]

	tau = size(p)[1]

	energy_stat = sum((p.^2)/2 + v.(r), dims = 1)/tau

	return energy_stat[1, :]

end
