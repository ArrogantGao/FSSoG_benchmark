using FastSpecSoG, Plots, CSV, DataFrames, ExTinyMD, LaTeXStrings
using Random
Random.seed!(123)

n_atoms = 1000
L = (20.0, 20.0, 20.0)

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]


poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

uspara_array = [USeriesPara(i...) for i in 3:6]
for i in 1:4
    uspara = uspara_array[i]
    M_mid = length(uspara.sw)
    E_kspace=[]

    #E_direct_ks = long_energy_us_k(qs, poses, 50, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)
    E_direct_true = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)

    for kk in 1:20
        E_temp = long_energy_us_k(qs, poses, kk, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)
        error = abs(E_direct_true-E_temp)
        error_rel = error/abs(E_direct_true)
        push!(E_kspace, error)

        df = DataFrame(uspara = i+2, energy_exact = E_direct_true, abs_error = error, relative_error = error_rel)
        CSV.write("data_for_manu/Acc_T1_Kspace.csv", df, append = true)
    end
end
