using FastSpecSoG, StatProfilerHTML, BenchmarkTools
using Plots
using Random
Random.seed!(1234)

n_atoms = 36

L = (100.0, 100.0, 100.0)

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

# using the FFCT

N_reals = [(32, 32, 32), (64, 64, 64), (128, 128, 128)]
# w = Int.(N_real ./ 4)
w = (8, 8, 8)
β = 5.0 .* w
extra_pad_ratio = 2
cheb_order = 10
uspara = USeriesPara(5)

errors = []

M_mid_array = [1:2:20...]

E_direct_ks = [long_energy_us_k(qs, poses, 50, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid) for M_mid in M_mid_array]

for N_real in N_reals
    error_N = []
    @show N_real
    for i in 1:length(M_mid_array)
        M_mid = M_mid_array[i]
        @show M_mid

        gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
        E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

        error = abs(E_3DFFT - E_direct_ks[i])

        @show M_mid, error

        push!(error_N, error)
    end
    push!(errors, error_N)
end

plot(xlabel="M_mid", ylabel="Error", title="Error vs M_mid, w = (8, 8, 8)", dpi = 500, ylim = [-12, 0])
plot!(M_mid_array, log10.(errors[1]), label="N_real = $(N_reals[1])", marker = :circle)
plot!(M_mid_array, log10.(errors[2]), label="N_real = $(N_reals[2])", marker = :square)
plot!(M_mid_array, log10.(errors[3]), label="N_real = $(N_reals[3])", marker = :diamond)
savefig("figs/FFT3D_accuracy_M_mid.png")