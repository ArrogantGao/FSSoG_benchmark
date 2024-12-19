using FastSpecSoG, CSV, DataFrames
using Random
Random.seed!(123)

@show Threads.nthreads(), nprocs()

n_atoms = 1000
L = (20.0, 20.0, 20.0)

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
uspara_array = [USeriesPara(i...) for i in 1:6]
uspara = uspara_array[6]
M_mid_array = [1, 3, 5, 7, 9]
N_xy_array = [1:1:16...]

data_file = "data/Acc_T4_cube_Nxy.csv"
#E_direct_total = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, length(uspara.sw))
E_direct_total = 293.04310539443907

E_ks = [166.05167919834435, 104.58856911741812, 60.897155401597864, 31.39401300686976, 13.414095784168726]

for i in 1:length(M_mid_array)
    M_mid = M_mid_array[i]
    # E_direct_k = long_energy_us_k(qs, poses, 1e-16, L, uspara, M_mid + 1, length(uspara.sw))
    E_direct_k = E_ks[i]

    for N_xy in N_xy_array
        N_grid = (N_xy, N_xy, 64)

        k_x, k_y, r_z, H_r, H_c, phase_x, phase_y = long_paras_gen(L, N_grid)
        cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, 64)
        E_FFCT_loop_k = energy_long_cheb_k(qs, poses, L, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, cheb_mat)

        abs_error = abs(E_FFCT_loop_k - E_direct_k)
        relative_error = abs(E_FFCT_loop_k - E_direct_k) / abs(E_direct_total)

        @show N_grid, M_mid, abs_error, relative_error

        df = DataFrame(N_xy = N_xy, M_mid = M_mid, E_exact = E_direct_k, E_FFCT = E_FFCT_loop_k, abs_error = abs_error, relative_error = relative_error)

        CSV.write(data_file, df, append=true)
    end
end