using FastSpecSoG, Plots, CSV, DataFrames, ExTinyMD, LaTeXStrings
using Random
Random.seed!(123)

n_atoms = 100
L = (100.0, 100.0, 100.0)

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

uspara_array = [USeriesPara(i...) for i in 1:6]
for i in 1:6
    uspara = uspara_array[i]
    M_mid = length(uspara.sw)
    E_kspace=[]

    #E_direct_ks = long_energy_us_k(qs, poses, 50, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)
    E_direct_true = long_energy_us_k(qs, poses, 100, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)

    for kk in 1:100
        E_temp = long_energy_us_k(qs, poses, kk, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)
        error = abs(E_direct_true-E_temp)
        error_rel = error/abs(E_direct_true)
        push!(E_kspace, error)

        df = DataFrame(uspara = i, energy_exact = E_direct_true, abs_error = error, relative_error = error_rel)
        CSV.write("data/Acc_T1_1_Kspace.csv", df, append = true)
    end
end
