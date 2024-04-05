using FastSpecSoG, CSV, DataFrames
using Random
Random.seed!(123)

n_atoms = 1000
L = (100.0, 100.0, 1.0)

qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
uspara = USeriesPara(6)
M_total = length(uspara.sw)

# E_direct_k = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_total)
# E_direct_0 = long_energy_us_0(qs, poses, L, uspara, 1, M_total)
# E_direct_total = E_direct_k + E_direct_0
# @show E_direct_k, E_direct_0, E_direct_total

E_direct_k, E_direct_0, E_direct_total = [279.1547597171813, 0.03041146481865265, 279.185171182]

N_real = (256, 256)
R_zs = [1:16...]
w = (16, 16)

errors = []
for R_z in R_zs
    β = 7.5 .* w
    cheb_order = 32
    Taylor_Q = 16

    gridinfo, pad_grids, cheb_coefs, scalefactors, H_r, H_c, cheb_value, r_z = thin_paras_gen(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)

    E_thin_k = energy_long_thin_k(qs, poses, L, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value)

    error_rel = abs(E_thin_k - E_direct_k) / abs(E_direct_total)

    df = DataFrame(R_z = R_z, E_exact = E_direct_k, E_long = E_thin_k, error_rel = error_rel)
    CSV.write("data/Acc_T5_thin_Rz.csv", df, append = true)
    push!(errors, error_rel)
    @show R_z, error_rel
end