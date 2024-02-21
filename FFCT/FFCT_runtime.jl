using FastSpecSoG, StatProfilerHTML, BenchmarkTools, Plots

n_atoms_array = [10^i for i in 1:5]
t_loop = []
t_einsum = []

for n_atoms in n_atoms_array

    L = (100.0, 100.0, 100.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    N_grid = (8, 8, 10)
    uspara = USeriesPara(2)
    soepara = SoePara4()
    M_mid = 8

    k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid, USeriesPara(2), M_mid, n_atoms)

    time_einsum = @belapsed energy_long_einsum($qs, $poses, $L, $M_mid, $k_x, $k_y, $k_mat, $r_z, $phase_x, $phase_y, $phase_xs, $phase_ys, $phase_xys, $temp_ijlk, $temp_ijl, $size_dict, $z_coef, $exp_coef, $z, $sort_z, $us_mat, $b_l, $b_u, $rhs, $sol, $ivsm, $H_r, $H_c, $H_s, $uspara, $soepara)
    push!(t_einsum, time_einsum)

    time_loop = @belapsed energy_long_loop($qs, $poses, $L, $N_grid, $M_mid, $k_x, $k_y, $r_z, $phase_x, $phase_y, $z, $sort_z, $us_mat, $b_l, $b_u, $rhs, $sol, $ivsm, $H_r, $H_c, $H_s, $uspara, $soepara)
    push!(t_loop, time_loop)

    @show n_atoms, time_einsum, time_loop
end


plot(xlabel = "log10(n_atoms)", ylabel = "log10(time)", title = "FFCT runtime comparison", legend = :topleft)
plot!([1:5...], log10.(t_einsum), label = "einsum", marker = :circle)
plot!([1:5...], log10.(t_loop), label = "loop", marker = :circle)
savefig("./figs/FFCT_runtime.png")