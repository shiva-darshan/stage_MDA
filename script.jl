using BenchmarkTools
include("./atom_chains.jl")
a, f_0, alpha  = 1., 1., 1.
theta = 1.
gamma_l, gamma_r = 2., 2.
lambda = 1.
T_l, T_r = 1., 1.

T, dt = 200_000., 0.005
skip = 10
N = floor(Int64, T/dt)
burn = 1.

n = 100
recenter = false
rotar = false

sampling_interval = 1
ind_traj = 1:n

ind_per = 1:n
period = 100.

q_bins = -4:0.05:4
p_bins = -4:0.05:4
r_bins = -4:0.05:4

return_counts = false



F = get_forcing_cos(f_0, 1/theta, n)
tension = 1.
#F(t) = 0 

q_0 = zeros(n)
p_0 = zeros(n)

U, grad_U!, v, v_p, v_s = get_pot_fput(alpha, a, 1.)



@btime trans_coef, avg_trans_coef = atom_chain_trans_coef(q_0, p_0, grad_U!, v_p, F, T_l, T_r, gamma_l, 
	gamma_r, tension, lambda, T, N, burn, skip; rotar, recenter, flip_right = true)

@profview trans_coef, avg_trans_coef = atom_chain_trans_coef(q_0, p_0, grad_U!, v_p, F, T_l, T_r, gamma_l, 
gamma_r, tension, lambda, T, N, burn, skip; rotar, recenter, flip_right = true)


# @btime atom_chain(q_0, p_0, grad_U!, v, v_p, v_s, F, T_l, T_r, gamma_l, gamma_r,
# 	tension, lambda, T, N, burn; pinning = alpha, recenter = false)

# @profview atom_chain(q_0, p_0, grad_U!, v, v_p, v_s, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, 
# 	T, N, burn; pinning = alpha, recenter = false)

# _, _, tp1, tp2, tension, _, _, _, _, _ = atom_chain(q_0, p_0, grad_U!, v, v_p, v_s, F, T_l, T_r, gamma_l, gamma_r,
# 	tension, lambda, T, N, burn; pinning = alpha, recenter = false)

# q_traj = nothing
# p_traj = nothing

# @btime atom_chain_traj(q_0, p_0, grad_U!, F, T_l, T_r, gamma_l, gamma_r, lambda, T, N, burn; 
# 	recenter = false)

# @profview atom_chain_traj(q_0, p_0, grad_U!, F, T_l, T_r, gamma_l, gamma_r, lambda, T, N, burn; 
# 	recenter = false)

# q_traj, p_traj = atom_chain_traj(q_0, p_0, grad_U!, F, T_l, T_r, gamma_l, gamma_r, lambda, T, N, 
# 	burn; recenter = false)

# @btime atom_chain_hist(q_0, p_0, grad_U!, F, T_l, T_r, gamma_l, gamma_r, lambda, T, N, burn, period, q_bins,
# 	p_bins, r_bins; recenter = false)

# @profview atom_chain_hist(q_0, p_0, grad_U!, F, T_l, T_r, gamma_l, gamma_r, lambda, T, N, burn, period, q_bins,
# 	p_bins, r_bins; recenter = false)

# @btime atom_chain_center(q_0, p_0, grad_U!, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, T, N, 
# 	burn, 10; recenter,  flip_right = true)

# @profview atom_chain_center(q_0, p_0, grad_U!, F, T_l, T_r, gamma_l, gamma_r, tension, lambda, T, N, 
# 	burn, 10; recenter,  flip_right = true)