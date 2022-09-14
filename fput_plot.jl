using Plots

a = -1.6667
b = .4444

plot3b = plot(x -> x^2/2 + a*x^3/3 + b*x^4/4, -2, 5, ylims = (-2, 10), legend = false,
    title = string("FPUT Potential with a = ", a, " and b = ", b))
savefig(plot3b, "FPUT_pot_min3.png")