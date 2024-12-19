using FastSpecSoG, CSV, DataFrames, LaTeXStrings, ChebParticleMesh, FFTW, Interpolations
using Random
Random.seed!(123)

CSV.write("data/Acc_T5_thin_w_diffkernels.csv", DataFrame(preset = Int[], w = Int[], gamma = Int[], energy_exact = Float64[], energy_k = Float64[], abs_error = Float64[], relative_error = Float64[]))

include("diff_kernels.jl")

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
R_z = 16
ws = [(i, i) for i in 1:16]

for w in ws
    β = 7.5 .* w
    alpha = π.*0.91 .* w 
    cheb_order = 32
    Taylor_Q = 16

    gridinfo, pad_grids, cheb_coefs, scalefactors, H_r, H_c, cheb_value, r_z = thin_paras_gen(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)
    E_thin_k_pkb = energy_long_thin_k(qs, poses, L, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value)

    gridinfo, pad_grids, cheb_coefs, scalefactors = thin_grids_gen_ES(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)
    E_thin_k_ES = energy_long_thin_k(qs, poses, L, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value)

    gridinfo, pad_grids, cheb_coefs, scalefactors = thin_grids_gen_Gaussian(N_real, R_z, w, alpha, L, cheb_order, uspara, Taylor_Q)
    E_thin_k_Gaussian = energy_long_thin_k(qs, poses, L, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value)


    error_rel_pkb = abs(E_thin_k_pkb - E_direct_k) / abs(E_direct_total)
    error_rel_ES = abs(E_thin_k_ES - E_direct_k) / abs(E_direct_total)
    error_rel_Gaussian = abs(E_thin_k_Gaussian - E_direct_k) / abs(E_direct_total)


    df = DataFrame(w = w[1], error_rel_pkb = error_rel_pkb, error_rel_ES = error_rel_ES, error_rel_Gaussian = error_rel_Gaussian)
    CSV.write("data/Acc_T5_thin_w_diffkernel.csv", df, append = true)
    @show w[1], error_rel_pkb, error_rel_ES, error_rel_Gaussian
end