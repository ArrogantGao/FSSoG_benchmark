using FastSpecSoG, CSV, DataFrames
using Random
Random.seed!(123)

CSV.write("data/Acc_T5_thin_Nxy.csv", DataFrame(preset = Int[], N_xy = Int[], gamma = Int[], energy_exact = Float64[], energy_k = Float64[], abs_error = Float64[], relative_error = Float64[]))

n_atoms = 1000

qs = [(-1.0)^i for i in 1:n_atoms]
uspara = USeriesPara(6)
M_total = length(uspara.sw)

# E_direct_k = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_total)
# E_direct_0 = long_energy_us_0(qs, poses, L, uspara, 1, M_total)
# E_direct_total = E_direct_k + E_direct_0
# @show E_direct_k, E_direct_0, E_direct_total

R_z = 16
cheb_order = 32
Taylor_Q = 16
Ns = Int.(ceil.([2^i for i in 4:0.5:10]))

gammas = [100, 500, 1000]
Es = [[279.15475971718234, 0.030411464763084064, 279.1851711819454], [295.955350846526, 0.0005254851891625051, 295.95587633171516], [258.7378007882928, 0.0009864338858588206, 258.7387872221787]]

for gamma in gammas
    L = (100.0, 100.0, 100.0 / gamma)
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
    i = findfirst(x -> x == gamma, gammas)
    E_direct_k, E_direct_0, E_direct_total = Es[i]
    @show E_direct_0, E_direct_k, E_direct_total
    for N in Ns
        N_real = (N, N)
        w = Int.(ceil.(N_real ./ 16))
        β = 7.5 .* w
        
        gridinfo, pad_grids, cheb_coefs, scalefactors, H_r, H_c, cheb_value, r_z = thin_paras_gen(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)

        E_thin_k = energy_long_thin_k(qs, poses, L, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value)

        error_rel = abs(E_thin_k - E_direct_k) / abs(E_direct_total)

        df = DataFrame(Nxy = N, gamma = gamma, E_exact = E_direct_k, E_long = E_thin_k, error_rel = error_rel)
        CSV.write("data/Acc_T5_thin_Nxy.csv", df, append = true)
        @show N, error_rel
    end
end