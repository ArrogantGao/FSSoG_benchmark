using FastSpecSoG, Plots, CSV, DataFrames, LaTeXStrings

using Random
Random.seed!(123)

CSV.write("data/Acc_T2_2_Extrapadding.csv", DataFrame(preset = Int[], extra_pad_ratio = Int[], M_mid = Int[], error = Float64[], error_rel = Float64[]))

n_atoms = 1000
L = (20.0, 20.0, 20.0)
N_real = (64, 64, 64)
w = (12, 12, 12)
β = 5.0 .* w
extra_pad_ratios = [0:30...]
cheb_order = 10
uspara = USeriesPara(6)
M_mids = [10 14 18 22]

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

errors = [[] for i in 1:length(M_mids)]
errors_rel = [[] for i in 1:length(M_mids)]

# E_direct_ks = [long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid) for M_mid in M_mids]
E_direct_ks = [224.75435272885954, 252.23674341893863, 268.4221529952294, 278.6463516342223]
@show E_direct_ks

#E_direct_ks_total = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, length(uspara.sw)) + long_energy_us_0(qs, poses, L, uspara, 1, length(uspara.sw)) 
E_direct_ks_total = 293.04310539443907

for i in 1:length(M_mids)
    for j in 1:length(extra_pad_ratios)
        extra_pad_ratio = extra_pad_ratios[j]
        M_mid = M_mids[i]

        gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
        E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

        error = abs(E_3DFFT - E_direct_ks[i])
        
        error_rel = error/ abs(E_direct_ks_total)

        @show extra_pad_ratio, M_mid, error, error_rel

        push!(errors[i], error)
        push!(errors_rel[i], error_rel)
        df = DataFrame(extra_pad_ratio = extra_pad_ratio, M_mid = M_mid, error = error, error_rel = error_rel)
        CSV.write("data/Acc_T2_2_Extrapadding.csv", df, append = true)
    end
end

