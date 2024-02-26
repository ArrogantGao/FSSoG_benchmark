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
        k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid, USeriesPara(2), M_mid, n_atoms)

        E_einsum = energy_long_einsum_k(qs, poses, L, M_mid, k_x, k_y, k_mat, r_z, phase_x, phase_y, phase_xs, phase_ys, phase_xys, temp_ijlk, temp_ijl, size_dict, z_coef, exp_coef, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, uspara)

        error = abs(E_einsum - E_direct_ks[M_mid])

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