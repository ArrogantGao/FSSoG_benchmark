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

uspara = USeriesPara(2)

errors = []

E_direct_ks = [long_energy_us_k(qs, poses, 30, L, uspara, M_mid + 1, length(uspara.sw)) for M_mid in 1:16]

for N_grid in [(4, 4, 4), (8, 8, 8), (16, 16, 16)]
    error_N = []
    @show N_grid
    for M_mid in 1:16
        k_x, k_y, r_z, H_r, H_c, phase_x, phase_y,  sort_z, z = long_paras_gen(L, N_grid, n_atoms)

        E_FFCT_loop_k = energy_long_loop_k(qs, poses, L, M_mid, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, uspara)

        error = abs(E_FFCT_loop_k - E_direct_ks[M_mid]) / abs(E_direct_ks[M_mid])

        @show M_mid, error

        push!(error_N, error)
    end
    push!(errors, error_N)
end

plot(xlabel = "M_mid", ylabel = "log10(error)", title = "FFCT relative error", dpi = 500)
plot!([1:16...], log10.(errors[1]), label = "N_grid = (4, 4, 4)", marker = :circle)
plot!([1:16...], log10.(errors[2]), label = "N_grid = (8, 8, 8)", marker = :square)
plot!([1:16...], log10.(errors[3]), label = "N_grid = (16, 16, 16)", marker = :diamond)
savefig("./figs/FFCT_accuracy.png")