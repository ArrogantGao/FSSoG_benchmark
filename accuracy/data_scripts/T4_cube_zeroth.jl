using FastSpecSoG, CSV, DataFrames
using Random
Random.seed!(123)

CSV.write("data/Acc_T4_cube_zeroth.csv", DataFrame(preset = Int[], Q_0 = Int[], M_mid = Int[], energy_exact = Float64[], energy_k = Float64[], abs_error = Float64[], relative_error = Float64[]))

@show Threads.nthreads()

n_atoms = 1000
L = (20.0, 20.0, 20.0)

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
uspara_array = [USeriesPara(i...) for i in 1:6]
uspara = uspara_array[6]
M_mid_array = [1, 3, 5, 7, 9]
Q0_array = [1:4:64...]

# E_naives = [long_energy_us_0(qs, poses, L, uspara, M_mid + 1, length(uspara.sw)) for M_mid in M_mid_array]
E_naives = [87.05908690289544, 83.96166551720789, 79.49334643724278, 73.24642894038898, 65.06713063525174]
@show E_naives

for m in 1:length(M_mid_array)
    M_mid = M_mid_array[m]
    E_naive = E_naives[m]
    for Q_0 in Q0_array
        r_z, grids, chebcoefs = zero_paras_gen(L[3], Q_0)
        E_FGT = zeroth_order(qs, poses, L, uspara, M_mid, r_z, chebcoefs, grids)

        abs_error = abs(E_FGT - E_naive)
        relative_error = abs(E_FGT - E_naive) / abs(E_naive)

        @show Q_0, M_mid, E_FGT, E_naive, abs_error, relative_error

        df = DataFrame(Q_0 = Q_0, M_mid = M_mid, energy_exact = E_naive, energt_fgt = E_FGT, error_rel = relative_error)
        CSV.write("data/Acc_T4_cube_zeroth.csv", df, append=true)
    end
end