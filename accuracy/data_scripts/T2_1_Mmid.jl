using FastSpecSoG, CSV, DataFrames

using Random
Random.seed!(123)

CSV.write("data/Acc_T2_1_Mmid.csv", DataFrame(preset = Int[], M_mid = Int[], energy_exact = Float64[], energy_k = Float64[], abs_error = Float64[], relative_error = Float64[]))

n_atoms = 1000
L = (20.0, 20.0, 20.0)
N_real = (64, 64, 64)
w = (12, 12, 12)
β = 5.0 .* w
extra_pad_ratios = [2, 3, 4, 8]
cheb_order = 10
uspara = USeriesPara(6)
M_mids = [3:1:30...]

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

errors = [[] for i in 1:length(extra_pad_ratios)]
errors_rel = [[] for i in 1:length(extra_pad_ratios)]

# E_direct_ks = [long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid) for M_mid in M_mids]
E_direct_ks = [104.49287100862274, 130.34667369287365, 152.65260380440816, 171.87925994366861, 188.40266369599016, 202.53490948359013, 214.56187922382824, 224.75435272885954, 233.3684525923324, 240.65476933424884, 246.86489942455765, 252.23674341893863, 256.96854296656636, 261.19952788257365, 265.0061721212353, 268.4221529952294, 271.46831090953947, 274.16854444332273, 276.55116429233897, 278.6463516342223, 280.48410333677606, 282.09300048009175, 283.49955130466554, 284.72790274822864, 285.7997791195156, 286.7345544147295, 287.54939807167875, 288.25945625751814]
@show E_direct_ks

#E_direct_ks_total = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, length(uspara.sw)) + long_energy_us_0(qs, poses, L, uspara, 1, length(uspara.sw)) 
E_direct_ks_total = 293.04310539443907

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

