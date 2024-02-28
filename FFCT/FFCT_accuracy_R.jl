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

uspara = USeriesPara(5)
@show length(uspara.sw)

errors = []

M_mid_array = [1, 3, 6, 12, 24, 48]

E_direct_ks = [long_energy_us_k(qs, poses, 30, L, uspara, M_mid + 1, length(uspara.sw)) for M_mid in M_mid_array]

for N_grid in [(16, 16, j) for j in 2:60]
    error_N = []
    @show N_grid
    for i in 1:length(M_mid_array)
        M_mid = M_mid_array[i]
        @show M_mid

        k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid, uspara, M_mid, n_atoms)

        E_einsum = energy_long_einsum_k(qs, poses, L, M_mid, k_x, k_y, k_mat, r_z, phase_x, phase_y, phase_xs, phase_ys, phase_xys, temp_ijlk, temp_ijl, size_dict, z_coef, exp_coef, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, uspara)

        error = abs(E_einsum - E_direct_ks[i]) / abs(E_direct_ks[i])

        @show M_mid, error

        push!(error_N, error)
    end
    push!(errors, error_N)
end