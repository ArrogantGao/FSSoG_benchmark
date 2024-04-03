
using FastSpecSoG, CSV, DataFrames
using Plots
using Random
Random.seed!(123)

n_atoms = 1000

L = (20.0, 20.0, 20.0)

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]

poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

# using the FFCT

uspara = USeriesPara(6)

errors = []

M_mid_array = [1:20...]
N_z_array = [3:3:60...]


data_file = "data_for_manu/Acc_T4_2_Cheb.csv"
#E_direct_total = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, length(uspara.sw))
E_direct_total = 293.04310539443907

for i in 1:length(M_mid_array)
    M_mid = M_mid_array[i]
    E_direct_k = long_energy_us_k(qs, poses, 1e-16, L, uspara, M_mid + 1, length(uspara.sw))

    for N_z in N_z_array
        N_grid = (16, 16, N_z)    

        k_x, k_y, r_z, H_r, H_c, phase_x, phase_y,  sort_z, z = long_paras_gen(L, N_grid, n_atoms)

        E_FFCT_loop_k = energy_long_loop_k(qs, poses, L, M_mid, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, uspara)

        abs_error = abs(E_FFCT_loop_k - E_direct_k)
        relative_error = abs(E_FFCT_loop_k - E_direct_k) / abs(E_direct_total)

        @show N_grid, M_mid, abs_error, relative_error

        df = DataFrame(N_z = N_z, M_mid = M_mid, E_exact = E_direct_k, E_FFCT = E_FFCT_loop_k, abs_error = abs_error, relative_error = relative_error)

        CSV.write(data_file, df, append=true)
    end
end


