using FastSpecSoG, CSV, DataFrames
using Plots
using Random
Random.seed!(1234)

n_atoms = 100

L = (100.0, 100.0, 100.0)

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

# using the FFCT

uspara = USeriesPara(2)

errors = []

M_mid_array = [1:2:9...]
N_xy_array = [2, 4, 8, 16, 32, 64]


data_file = "data/accuracy_long_Nxy.csv"

for i in 1:length(M_mid_array)
    M_mid = M_mid_array[i]
    E_direct_k = long_energy_us_k(qs, poses, 30, L, uspara, M_mid + 1, length(uspara.sw))

    for N_xy in N_xy_array
        N_grid = (N_xy, N_xy, 32)    

        k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid, uspara, M_mid, n_atoms)

        E_FFCT_einsum_k = energy_long_einsum_k(qs, poses, L, M_mid, k_x, k_y, k_mat, r_z, phase_x, phase_y, phase_xs, phase_ys, phase_xys, temp_ijlk, temp_ijl, size_dict, z_coef, exp_coef, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, uspara)

        abs_error = abs(E_FFCT_einsum_k - E_direct_k)
        relative_error = abs(E_FFCT_einsum_k - E_direct_k) / abs(E_direct_k)

        @show N_grid, M_mid, abs_error, relative_error

        df = DataFrame(N_xy = N_xy, M_mid = M_mid, E_exact = E_direct_k, E_FFCT = E_FFCT_einsum_k, abs_error = abs_error, relative_error = relative_error)

        CSV.write(data_file, df, append=true)
    end
end
