using FastSpecSoG, Plots, CSV, DataFrames, ExTinyMD, LaTeXStrings
using Random
Random.seed!(123)

n_atoms = 1000
L = (20.0, 20.0, 20.0)

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
K_max_array = [3,6,9,12]
uspara_array = [USeriesPara(i...) for i in 1:6]
uspara = uspara_array[6]
M_mid_array = [1:1:15...]
M_total = length(uspara.sw)
E_kspace=[]
E_direct_total = long_energy_us_k(qs, poses, 20, L, uspara, 1, M_total) + long_energy_us_0(qs, poses, L, uspara, 1, M_total)
global E_direct_true = E_direct_total
@show E_direct_true

for j in 1:length(M_mid_array)
    M_mid = M_mid_array[j]
    #E_direct_ks = long_energy_us_k(qs, poses, 50, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)
    global E_direct_true -= long_energy_us_k(qs, poses, 20, L, uspara, M_mid, M_mid) + long_energy_us_0(qs, poses, L, uspara, M_mid, M_mid)
    @show E_direct_true

    for i in 1:length(K_max_array)
        kk = K_max_array[i]
        E_temp = long_energy_us_k(qs, poses, kk, L, uspara, M_mid+1, M_total) + long_energy_us_0(qs, poses, L, uspara, M_mid+1, M_total)
        error = abs(E_direct_true-E_temp)
        error_rel = error/abs(E_direct_total)
        push!(E_kspace, error)

        @show kk, M_mid, error, error_rel
        df = DataFrame(k_max = kk, energy_exact = E_direct_true, abs_error = error, relative_error = error_rel, M_mid = M_mid)
        CSV.write("data_for_manu/Acc_T4_1_wideKspace.csv", df, append = true)
    end
end

