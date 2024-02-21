using FastSpecSoG, StatProfilerHTML, BenchmarkTools

n_atoms = 1000

L = (100.0, 100.0, 100.0)

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

# using the 3DFFT
N_grid = (8, 8, 10)
uspara = USeriesPara(2)
soepara = SoePara4()
M_mid = 8

k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid, USeriesPara(2), M_mid, n_atoms)

energy_long_einsum(qs, poses, L, M_mid, k_x, k_y, k_mat, r_z, phase_x, phase_y, phase_xs, phase_ys, phase_xys, temp_ijlk, temp_ijl, size_dict, z_coef, exp_coef, z, sort_z, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, uspara, soepara)

@profilehtml for _=1:1000 energy_long_einsum(qs, poses, L, M_mid, k_x, k_y, k_mat, r_z, phase_x, phase_y, phase_xs, phase_ys, phase_xys, temp_ijlk, temp_ijl, size_dict, z_coef, exp_coef, z, sort_z, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, uspara, soepara) end