using FastSpecSoG, Plots, CSV, DataFrames, ExTinyMD, LaTeXStrings
using Interpolations, FFTW
using ChebParticleMesh
using Random
Random.seed!(123)

CSV.write("data/Acc_T3_2_Nspec_PKB.csv", DataFrame(preset = Int[], N_real = Int[], energy_exact = Float64[], energy_k = Float64[], abs_error = Float64[], relative_error = Float64[]))

n_atoms = 1000
L = (20.0, 20.0, 20.0)

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

uspara = USeriesPara(6)
M_mid = 8

#E_direct_true = long_energy_us_k(qs, poses, 20, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)
#@show E_direct_true

#n_atoms0 = (256, 256, 256)
#ratios = [2.0 ^i for i in 0.0:0.5:2.0]
N_real_array = [(i,i,i) for i in 12:4:64]
w = (12, 12, 12) 
β = 5.0 .* w       
alpha = π.*0.91 .* w 
@show N_real_array



function Tdk(kx::T, ky::T, kz::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    k2 = kx^2 + ky^2 + kz^2

    Tdk_value = zero(T)
    for i in 1:M_mid
        sl, wl = uspara.sw[i]
        Tdk_value += wl * sl^3 * exp(- sl^2 * k2 / 4)
    end

    return sqrt(π)/4 * Tdk_value
end


function mid_paras_gen_diffwindow(N_real::NTuple{3, Int}, w::NTuple{3, Int}, β::NTuple{3, T}, L::NTuple{3, T}, extra_pad_ratio::Int, cheb_order::Int, uspara::USeriesPara{T}, M_mid::Int, alpha::NTuple{3, T}) where{T}
    periodicity = (true, true, false)
    extra_pad = extra_pad_ratio .* w
    @assert M_mid <= length(uspara.sw)
    gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
    gridbox = GridBox(gridinfo)
    fw = [x -> Wkb(x, (w[i] + 0.5) * gridinfo.h[i], β[i]) for i in 1:3]
    F_fw = [x -> FWkb(x, (w[i] + 0.5) * gridinfo.h[i], β[i]) for i in 1:3]
    #H = [(w[i] + 0.5) * gridinfo.h[i] for i in 1:3]
    #fw = [x -> abs(x) < H[i] ? exp(-alpha[i] * (x/H[i])^2) : 0.0 for i in 1:3]
    #F_fw = [x -> sqrt(π/alpha[i]) * H[i] * exp(-x^2 * H[i]^2 / (4 * alpha[i]) )  for i in 1:3]

    cheb_coefs = tuple([ChebCoef(fw[i], gridinfo.h[i], w[i], cheb_order) for i in 1:3]...)
    func_scale = (kx, ky, kz) -> (F_fw[1](kx) * F_fw[2](ky) * F_fw[3](kz))^(-2) * Tdk(kx, ky, kz, uspara, M_mid)
    scalefactor = ScalingFactor(func_scale, gridinfo)
    return (gridinfo, gridbox, cheb_coefs, scalefactor)
end


extra_pad_ratio = 8
cheb_order = 10
# E_N_grid=[]
# N_real_true = (40, 40, 40)
# gridinfo_T, gridbox_T, cheb_coefs_T, scalefactor_T = mid_paras_gen_diffwindow(N_real_true, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid, alpha)
# E_direct_true = energy_mid(qs, poses, gridinfo_T, gridbox_T, cheb_coefs_T, scalefactor_T)

# E_direct_true = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)
E_direct_true =  202.53490948359013
@show E_direct_true


for N_real in N_real_array
    if mod(N_real[1], 2) == 1
        N_real = N_real .+ (1, 1, 1)
    end
        gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen_diffwindow(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid, alpha)
        E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

        error = abs(E_3DFFT - E_direct_true)
        error_rel = error / abs(E_direct_true)
      
        @show N_real, error, error_rel

        df = DataFrame(Window_type = 1, energy_exact = E_direct_true, abs_error = error, relative_error = error_rel, N_real = N_real[1])
        #CSV.write("data/Acc_T1_2_WinFunc_Gaussian.csv", df, append = true)
        CSV.write("data/Acc_T3_2_Nspec_PKB.csv", df, append = true)
end
