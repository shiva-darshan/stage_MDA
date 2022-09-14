"""
    Functions for compute various quanties from the analytic experessions
    in the article KLO 2022

"""

function green_func(x, alpha, theta, l, gamma, n)
    summand = 0. *im
    
    for j = 0:n
        num = (-1)^j * cos((pi*j*(2*x+1))/(2*(n+1))) * cos(pi*j/(2*(n+1)))
        denom = 4*sin(pi*j/(2*(n+1)))^2 + alpha - (2*pi*l/theta)^2 + (4*pi*gamma*l/theta)*im
        summand += (2 - (j==0))/(n+1) * (num/denom)
    end
    
    return summand
end


function q_tilde(l, x, alpha, theta, gamma, n, c)
    green = green_func(x, alpha, theta, l, gamma, n)
    
    return abs(l) == 1 ? 0.5*c*green/sqrt(n+1) : 0. 
end

function sum_q_tilsq(l, alpha, theta, gamma, n, c)
    summand(j) = (2 - (j == 0)) * (cos(j*pi/(2*(n+1)))^2)/(
        (4*sin(j*pi/(2*(n+1)))^2 + alpha - (2*pi*l/theta)^2)^2 + (4*gamma*pi*l/theta)^2)

    return sum(summand.(0:n))/(n+1)^2 *(c/2)^2
end

function J_n(alpha, theta, gamma, n, c)
    
    return -4*gamma*(2*pi/theta)^2 * sum_q_tilsq(1, alpha, theta, gamma, n, c)
end

function p_tilde(l, x, alpha, theta, gamma, n, c)
    return q_tilde(l, x, alpha, theta, gamma, n, c)*((2*pi*im*l)/theta)
end

function cycle_shift(x, k)
    return [x[(k+1):end]; x[1:k]]        
end

D(alpha) = 2/(2 + alpha + sqrt(alpha^2 + 4*alpha))


"""
Some stuff for computing local equlibirum stats
"""


function LTE_energy_stretch(temps, a,b)
    d = length(temps)
    energies = zeros(d)
    r = zeros(d-1)

    v(r) = r^2/2 + a*r^3/3 + b*r^4/4

    energies[1] = temps[1]/2

    for i = 2:d
        T_eff = temps[i]
        M = 2.5
        M = max(M/T_eff, M * T_eff)
        xs = LinRange(-10, 10, 100001)
        Z = sum((exp.(-v.(xs[1:end-1])/T_eff) + exp.(-v.(xs[2:end])/T_eff))/2)*(xs[2] - xs[1])
        norm_pdf(x) = exp(-(x -1)^2/2)/sqrt(2 * pi)

        #density of FPUT bond diameter probability
        f(x) = exp(-v(x)/T_eff)/Z


        energies[i] = sum((v.(xs[1:end-1]).*f.(xs[1:end-1]) 
            .+ v.(xs[2:end]).*f.(xs[2:end]))/2)*(xs[2] - xs[1]) + T_eff/2
        r[i-1] = sum((xs[1:end-1].*f.(xs[1:end-1]) .+ xs[2:end].*f.(xs[2:end]))/2)*(xs[2] - xs[1])
    end

    return energies, r
end

function LTE_energy_stretch_cos(temps)
    d = length(temps)
    energies = zeros(d)
    r = zeros(d-1)

    v(r) = 1 - cos(r)

    energies[1] = temps[1]/2

    for i = 2:d
        T_eff = temps[i]
        M = 2.5
        M = max(M/T_eff, M * T_eff)
        xs = LinRange(-pi, pi, 10001)
        Z = sum((exp.(-v.(xs[1:end-1])/T_eff) + exp.(-v.(xs[2:end])/T_eff))/2)*(xs[2] - xs[1])
        norm_pdf(x) = exp(-(x -1)^2/2)/sqrt(2 * pi)

        #density of cos bond diameter probability
        f(x) = exp(-v(x)/T_eff)/Z


        energies[i] = sum((v.(xs[1:end-1]).*f.(xs[1:end-1]) 
            .+ v.(xs[2:end]).*f.(xs[2:end]))/2)*(xs[2] - xs[1]) + T_eff/2
        r[i-1] = sum((xs[1:end-1].*f.(xs[1:end-1]) .+ xs[2:end].*f.(xs[2:end]))/2)*(xs[2] - xs[1])
    end

    return energies, r
end