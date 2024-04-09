using FastSpecSoG, CSV, DataFrames
using Random
Random.seed!(123)

n_atoms = 1000

qs = [(-1.0)^i for i in 1:n_atoms]
uspara = USeriesPara(6)
M_total = length(uspara.sw)

# E_direct_k = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_total)
# E_direct_0 = long_energy_us_0(qs, poses, L, uspara, 1, M_total)
# E_direct_total = E_direct_k + E_direct_0
# @show E_direct_k, E_direct_0, E_direct_total

N_real = (256, 256)
R_zs = [1:16...]
w = (16, 16)
gammas = [100, 500, 1000]

M_mid = 0

E_direct_0s = [0.03041146499605954, 0.0005254852391897013, 0.0009864339577972327]

for i in 1:3
    gamma = gammas[i]
    L = (100.0, 100.0, 100.0 / gamma)
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
    # E_direct_0 = long_energy_us_0(qs, poses, L, uspara, 1, M_total)
    E_direct_0 = E_direct_0s[i]
    @show E_direct_0
    for R_z in R_zs
        β = 7.5 .* w
        cheb_order = 32
        Taylor_Q = 16

        r_z, grids, chebcoefs = zero_paras_gen(L[3], R_z)
        E_FGT = zeroth_order(qs, poses, L, uspara, M_mid, r_z, chebcoefs, grids)

        error_rel = abs(E_FGT - E_direct_0) / abs(E_direct_0)

        df = DataFrame(R_z = R_z, gamma=gamma, E_exact = E_direct_0, E_fgt = E_FGT, error_rel = error_rel)
        CSV.write("data/Acc_T5_thin_zeroth.csv", df, append = true)
        @show gamma, R_z, error_rel
    end
end