using FastSpecSoG, Plots, CSV, DataFrames, LaTeXStrings

using Random
Random.seed!(1234)

n_atoms = 100
L = (100.0, 100.0, 100.0)
N_real = (256, 256, 256)
w = (16, 16, 16)
β = 5.0 .* w
extra_pad_ratios = [2, 3, 4, 8]
cheb_order = 10
uspara = USeriesPara(6)
M_mids = [1:50...]

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

errors = [[] for i in 1:length(extra_pad_ratios)]
errors_rel = [[] for i in 1:length(extra_pad_ratios)]

E_direct_ks = [long_energy_us_k(qs, poses, 100, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid) for M_mid in M_mids]

E_direct_ks_total = long_energy_us_k(qs, poses, 100, L, uspara, 1, length(uspara.sw)) + long_energy_us_0(qs, poses, L, uspara, 1, length(uspara.sw)) 

for i in 1:length(extra_pad_ratios)
    for j in 1:length(M_mids)
        extra_pad_ratio = extra_pad_ratios[i]
        M_mid = M_mids[j]

        gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
        E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

        error = abs(E_3DFFT - E_direct_ks[j])
        
        error_rel = error/ abs(E_direct_ks_total)

        @show extra_pad_ratio, M_mid, error, error_rel

        push!(errors[i], error)
        push!(errors_rel[i], error_rel)
    end
end

df = DataFrame(extra_pad_ratio = repeat(extra_pad_ratios, inner = length(M_mids)), M_mid = repeat(M_mids, outer = length(extra_pad_ratios)), error = vcat(errors...), error_rel = vcat(errors_rel...))
CSV.write("data/Acc_T2_1_Mmid.csv", df)

