using LinearAlgebra
using LinearAlgebra.BLAS
using Random
using CUDA
using FileIO
using Profile
using BenchmarkTools


T, dt = 20000., 0.005
N = floor(Int64, T/dt)
n = 100

q_n = zeros(n)
p_n = zeros(n)
x_n = [q_n; p_n]


function syr_cpu!(q_n, p_n, N, n)
	S_mat = zeros(2n, 2n)
	x_n = [q_n; p_n]
	for i = 1:N
		q_n .= randn(n)
		p_n .= randn(n)
		x_n[1:n] .= q_n
		x_n[n+1:end] .= p_n

		syr!('U', 1., x_n, S_mat)
	end

	return nothing
end

function syr_cuda!(q_n, p_n, N, n)
	S_mat = CuArray{Float64}(undef, (2n, 2n))
	x_n = CuArray([q_n; p_n])
	CUDA.allowscalar(true)
	for i = 1:N
		q_n .= randn(n)
		p_n .= randn(n)
		for j = 1:n
			x_n[j] = q_n[j]
			x_n[n+j] = p_n[j]
		end

		CUBLAS.syr!('U', 1., x_n, S_mat)
	end

	return nothing
end

function cpu_test(q_n, p_n, N, n)
	@btime syr_cpu!(q_n, p_n, N, n)
	flush(stdout)

	Profile.clear()
	@profile syr_cpu!(q_n, p_n, N, n)
	save("cpu_test.jlprof",  Profile.retrieve()...)

	return nothing
end

function cuda_test(q_n, p_n, N, n)
	@btime syr_cuda!(q_n, p_n, N, n)
	flush(stdout)

	Profile.clear()
	@profile syr_cuda!(q_n, p_n, N, n)
	save("cuda_test.jlprof",  Profile.retrieve()...)

	return nothing
end