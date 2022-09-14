using JLD
using Plots
using Statistics
include("../code/outils.jl")
include("../code/atom_chains.jl")


datadir = "../code/data/forced_fput/"
c = 10.0
alpha = 0.
omega = 1.
T_l = 0.1
gamma_l = 2.
lambda = 1.
tension = 0.
dt = 0.005
period = 0.
burn = 1.0e6

params(seed, d, a, b, T, burn) = string("_seed", seed, "_a",a, "_b", b, "_c", c, "_alpha", alpha, 
	"_omega", omega, "_Tl", T_l, "_gamma_l", gamma_l, "_lambda", lambda, "_tension", tension, "_d", 
	d, "_T", T, "_burn", burn, "_dt", dt, "_period", period, ".jld")


seed = 9104973537
a = 5.
b = 4.
T = 2.5e8
d = 1000

M_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
M_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
M_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
M_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
M_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
M_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")

M_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
M_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
M_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
M_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
M_mean_tau = (M_tau1 + M_tau2)/2
M_mean_r = (M_r1 + M_r2)/2
M_mean_temp = (M_temp1 + M_temp2)/2
M_mean_j = (M_j_1 + M_j_2)/2
M_mean_eng = (M_eng1 + M_eng2)/2

seed = 7408454679
d = 500

D_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
D_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
D_tau3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "tension_profile")
D_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
D_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
D_eng3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "e_profile")
D_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
D_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")
D_r3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "r_profile")

D_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
D_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
D_j_3 = mean(load(datadir*"cont2FPUT_fluxes"*params(seed, d, a, b, 2*T, 0.), "e_flux")[2:end])
D_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
D_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
D_temp3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "T_profile")
D_mean_tau = 0.25*D_tau1 + 0.25*D_tau2 + 0.5*D_tau3
D_mean_r = 0.25*D_r1 + 0.25*D_r2 + 0.5*D_r3
D_mean_temp = 0.25*D_temp1 + 0.25*D_temp2 + 0.5*D_temp3
D_mean_j = 0.25*D_j_1 + 0.25*D_j_2 + 0.5*D_j_3
D_mean_eng = 0.25*D_eng1 + 0.25*D_eng2 + 0.5*D_eng3


seed = 9080822348
T = 2.5e7
d = 100

C_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
C_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
C_tau3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "tension_profile")
C_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
C_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
C_eng3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "e_profile")
C_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
C_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")
C_r3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "r_profile")

C_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
C_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
C_j_3 = mean(load(datadir*"cont2FPUT_fluxes"*params(seed, d, a, b, 2*T, 0.), "e_flux")[2:end])
C_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
C_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
C_temp3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "T_profile")
C_mean_tau = 0.25*C_tau1 + 0.25*C_tau2 + 0.5*C_tau3
C_mean_r = 0.25*C_r1 + 0.25*C_r2 + 0.5*C_r3
C_mean_temp = 0.25*C_temp1 + 0.25*C_temp2 + 0.5*C_temp3
C_mean_j = 0.25*C_j_1 + 0.25*C_j_2 + 0.5*C_j_3
C_mean_eng = 0.25*C_eng1 + 0.25*C_eng2 + 0.5*C_eng3

p1 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300))
plot!(p1, C_mean_temp[1:end-1], -C_mean_j ./ (C_mean_temp[2:end] - C_mean_temp[1:end-1]), 
	label = "100 atoms")
p2 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300), 
    title = "Thermal Conductivity with Zero Tension")
plot!(p2, D_mean_temp[1:end-1], -D_mean_j ./ (D_mean_temp[2:end] - D_mean_temp[1:end-1]), 
	label = "500 atoms")
p3 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300))
plot!(p3, M_mean_temp[1:end-1], -M_mean_j ./ (M_mean_temp[2:end] - M_mean_temp[1:end-1]), 
	label = "1000 atoms")
plot1a = plot(p1, p2, p3, layout = (1, 3), size = (900, 300), margin=5Plots.mm)
savefig(plot1a, "therm_cond_min-1.png")

p1 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300))
plot!(p1, LinRange(0, 1, 100), C_mean_eng, label = "100 atoms")
p2 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300), 
    title = "Energy with Zero Tension")
plot!(p2, LinRange(0, 1, 500), D_mean_eng, label = "500 atoms")
p3 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300))
plot!(p3, LinRange(0, 1, 1000), M_mean_eng, label = "1000 atoms")
plot1b = plot(p1, p2, p3, layout = (1, 3), size = (900, 300), margin=5Plots.mm)
savefig(plot1b, "energy_min-1.png")

# plot1b = plot(x -> x^2/2 + a*x^3/3 + b*x^4/4, -2, 1, ylims = (-0.5, 2.5), legend = false,
#     title = string("FPUT Potential with a = ", a, " and b = ", b))
# savefig(plot1b, "FPUT_pot_min-1.png")



seed = 9104973537
a = -2.5
b = 1.
T = 2.5e8
d = 1000

M_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
M_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
M_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
M_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
M_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
M_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")

M_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
M_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
M_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
M_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
M_mean_tau = (M_tau1 + M_tau2)/2
M_mean_r = (M_r1 + M_r2)/2
M_mean_temp = (M_temp1 + M_temp2)/2
M_mean_j = (M_j_1 + M_j_2)/2
M_mean_eng = (M_eng1 + M_eng2)/2

seed = 7408454679
d = 500

D_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
D_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
D_tau3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "tension_profile")
D_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
D_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
D_eng3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "e_profile")
D_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
D_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")
D_r3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "r_profile")

D_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
D_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
D_j_3 = mean(load(datadir*"cont2FPUT_fluxes"*params(seed, d, a, b, 2*T, 0.), "e_flux")[2:end])
D_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
D_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
D_temp3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "T_profile")
D_mean_tau = 0.25*D_tau1 + 0.25*D_tau2 + 0.5*D_tau3
D_mean_r = 0.25*D_r1 + 0.25*D_r2 + 0.5*D_r3
D_mean_temp = 0.25*D_temp1 + 0.25*D_temp2 + 0.5*D_temp3
D_mean_j = 0.25*D_j_1 + 0.25*D_j_2 + 0.5*D_j_3
D_mean_eng = 0.25*D_eng1 + 0.25*D_eng2 + 0.5*D_eng3


seed = 9080822348
T = 2.5e7
d = 100

C_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
C_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
C_tau3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "tension_profile")
C_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
C_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
C_eng3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "e_profile")
C_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
C_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")
C_r3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "r_profile")

C_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
C_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
C_j_3 = mean(load(datadir*"cont2FPUT_fluxes"*params(seed, d, a, b, 2*T, 0.), "e_flux")[2:end])
C_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
C_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
C_temp3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "T_profile")
C_mean_tau = 0.25*C_tau1 + 0.25*C_tau2 + 0.5*C_tau3
C_mean_r = 0.25*C_r1 + 0.25*C_r2 + 0.5*C_r3
C_mean_temp = 0.25*C_temp1 + 0.25*C_temp2 + 0.5*C_temp3
C_mean_j = 0.25*C_j_1 + 0.25*C_j_2 + 0.5*C_j_3
C_mean_eng = 0.25*C_eng1 + 0.25*C_eng2 + 0.5*C_eng3

p1 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300))
plot!(p1, C_mean_temp[1:end-1], -C_mean_j ./ (C_mean_temp[2:end] - C_mean_temp[1:end-1]), 
	label = "100 atoms")
p2 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300), 
    title = "Thermal Conductivity with Zero Tension")
plot!(p2, D_mean_temp[1:end-1], -D_mean_j ./ (D_mean_temp[2:end] - D_mean_temp[1:end-1]), 
	label = "500 atoms")
p3 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300))
plot!(p3, M_mean_temp[1:end-1], -M_mean_j ./ (M_mean_temp[2:end] - M_mean_temp[1:end-1]), 
	label = "1000 atoms")
plot1a = plot(p1, p2, p3, layout = (1, 3), size = (900, 300), margin=5Plots.mm)
savefig(plot1a, "therm_cond_min2.png")

p1 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300))
plot!(p1, LinRange(0, 1, 100)[2:end], C_mean_eng[2:end], label = "100 atoms")
p2 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300), 
    title = "Energy with Zero Tension")
plot!(p2, LinRange(0, 1, 500)[2:end], D_mean_eng[2:end], label = "500 atoms")
p3 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300))
plot!(p3, LinRange(0, 1, 1000)[2:end], M_mean_eng[2:end], label = "1000 atoms")
plot1b = plot(p1, p2, p3, layout = (1, 3), size = (900, 300), margin=5Plots.mm)
savefig(plot1b, "energy_min2.png")

# plot2b = plot(x -> x^2/2 + a*x^3/3 + b*x^4/4, -2, 4, ylims = (-1, 10), legend = false,
#     title = string("FPUT Potential with a = ", a, " and b = ", b))
# savefig(plot2b, "FPUT_pot_min2.png")



seed = 9104973537
a = -1.6667
b = .4444
T = 2.5e8
d = 1000

M_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
M_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
M_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
M_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
M_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
M_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")

M_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
M_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
M_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
M_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
M_mean_tau = (M_tau1 + M_tau2)/2
M_mean_r = (M_r1 + M_r2)/2
M_mean_temp = (M_temp1 + M_temp2)/2
M_mean_j = (M_j_1 + M_j_2)/2
M_mean_eng = (M_eng1 + M_eng2)/2

r_p = plot(legend = false)
plot!(r_p, M_mean_r)
savefig(r_p, "M_stretch_profile.png")

seed = 7408454679
d = 500

D_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
D_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
D_tau3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "tension_profile")
D_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
D_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
D_eng3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "e_profile")
D_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
D_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")
D_r3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "r_profile")

D_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
D_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
D_j_3 = mean(load(datadir*"cont2FPUT_fluxes"*params(seed, d, a, b, 2*T, 0.), "e_flux")[2:end])
D_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
D_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
D_temp3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "T_profile")
D_mean_tau = 0.25*D_tau1 + 0.25*D_tau2 + 0.5*D_tau3
D_mean_r = 0.25*D_r1 + 0.25*D_r2 + 0.5*D_r3
D_mean_temp = 0.25*D_temp1 + 0.25*D_temp2 + 0.5*D_temp3
D_mean_j = 0.25*D_j_1 + 0.25*D_j_2 + 0.5*D_j_3
D_mean_eng = 0.25*D_eng1 + 0.25*D_eng2 + 0.5*D_eng3

ten_p = plot(title = "Tension Profile", legend = :outertopright, size = (800, 400))
plot!(ten_p, D_tau1, label = string("Time = ", T))
plot!(ten_p, (D_tau1 + D_tau2)/2, label = string("Time = ", 2*T))
plot!(ten_p, D_mean_tau, label = string("Time = ", 4*T))
savefig(ten_p, "D_tension_profile.png")

r_p = plot(legend = false)
plot!(r_p, D_mean_r)
savefig(r_p, "D_stretch_profile.png")

seed = 9080822348
T = 2.5e7
d = 100

C_tau1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "tension_profile")
C_tau2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "tension_profile")
C_tau3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "tension_profile")
C_eng1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "e_profile")
C_eng2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "e_profile")
C_eng3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "e_profile")
C_r1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "r_profile")
C_r2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "r_profile")
C_r3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "r_profile")

C_j_1 = mean(load(datadir*"FPUT_fluxes"*params(seed, d, a, b, T, burn), "e_flux")[2:end])
C_j_2 = mean(load(datadir*"cont1FPUT_fluxes"*params(seed, d, a, b, T, 0.), "e_flux")[2:end])
C_j_3 = mean(load(datadir*"cont2FPUT_fluxes"*params(seed, d, a, b, 2*T, 0.), "e_flux")[2:end])
C_temp1 = load(datadir*"FPUT_profiles"*params(seed, d, a, b, T, burn), "T_profile")
C_temp2 = load(datadir*"cont1FPUT_profiles"*params(seed, d, a, b, T, 0.), "T_profile")
C_temp3 = load(datadir*"cont2FPUT_profiles"*params(seed, d, a, b, 2*T, 0.), "T_profile")
C_mean_tau = 0.25*C_tau1 + 0.25*C_tau2 + 0.5*C_tau3
C_mean_r = 0.25*C_r1 + 0.25*C_r2 + 0.5*C_r3
C_mean_temp = 0.25*C_temp1 + 0.25*C_temp2 + 0.5*C_temp3
C_mean_j = 0.25*C_j_1 + 0.25*C_j_2 + 0.5*C_j_3
C_mean_eng = 0.25*C_eng1 + 0.25*C_eng2 + 0.5*C_eng3

p1 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300))
plot!(p1, C_mean_temp[1:end-1], -C_mean_j ./ (C_mean_temp[2:end] - C_mean_temp[1:end-1]), 
	label = "100 atoms")
p2 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300), 
    title = "Energy with Zero Tension")
plot!(p2, D_mean_temp[1:end-1], -D_mean_j ./ (D_mean_temp[2:end] - D_mean_temp[1:end-1]), 
	label = "500 atoms")
p3 = plot(legend = :bottomright, xlabel = "Temperature", size = (300, 300))
plot!(p3, M_mean_temp[1:end-1], -M_mean_j ./ (M_mean_temp[2:end] - M_mean_temp[1:end-1]), 
	label = "1000 atoms")
plot1a = plot(p1, p2, p3, layout = (1, 3), size = (900, 300), margin=5Plots.mm)
savefig(plot1a, "therm_cond_min3.png")

p1 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300))
plot!(p1, LinRange(0, 1, 100)[2:end], C_mean_eng[2:end], label = "100 atoms")
p2 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300), 
    title = "Energy with Zero Tension")
plot!(p2, LinRange(0, 1, 500)[2:end], D_mean_eng[2:end], label = "500 atoms")
p3 = plot(legend = :bottomright, xlabel = "sites", size = (300, 300))
plot!(p3, LinRange(0, 1, 1000)[2:end], M_mean_eng[2:end], label = "1000 atoms")
plot1b = plot(p1, p2, p3, layout = (1, 3), size = (900, 300), margin=5Plots.mm)
savefig(plot1b, "energy_min3.png")

# plot3b = plot(x -> x^2/2 + a*x^3/3 + b*x^4/4, -2, 5, ylims = (-2, 10), legend = false,
#     title = string("FPUT Potential with a = ", a, " and b = ", b))
# savefig(plot3b, "FPUT_pot_min3.png")

plot4 = plot(ylims = (-2, 8), legend = :outertopright, title = "FPUT Potentials", size = (800, 400))
a, b = -2.5, 1.
plot!(plot4, x -> x^2/2 + a*x^3/3 + b*x^4/4, -2, 5, label = string("a = ", a, ", b = ", b))
a, b = -1.6667, 0.4444
plot!(plot4, x -> x^2/2 + a*x^3/3 + b*x^4/4, -2, 5, label = string("a = ", a, ", b = ", b))
a, b = 5, 4
plot!(plot4, x -> x^2/2 + a*x^3/3 + b*x^4/4, -2, 5, label = string("a = ", a, ", b = ", b))
savefig(plot4, "FPUT_pots.png")