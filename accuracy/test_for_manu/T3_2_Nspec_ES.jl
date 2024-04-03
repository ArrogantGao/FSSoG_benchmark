using FastSpecSoG, Plots, CSV, DataFrames, ExTinyMD, LaTeXStrings
using Interpolations, FFTW
using ChebParticleMesh
using Random
Random.seed!(123)

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
    #fw = [x -> Wkb(x, (w[i] + 0.5) * gridinfo.h[i], β[i]) for i in 1:3]
    #F_fw = [x -> FWkb(x, (w[i] + 0.5) * gridinfo.h[i], β[i]) for i in 1:3]
    H = [(w[i] + 0.5) * gridinfo.h[i] for i in 1:3]
    #fw = [x -> exp(-alpha[i] * (x/H[i])^2) for i in 1:3]
    #F_fw = [x -> sqrt(π/alpha[i]) * H[i] * exp(-x^2 * H[i]^2 / (4 * alpha[i]) )  for i in 1:3]
    fw = [x -> abs(x) < H[i] ? exp(β[i] * sqrt(1 - (x / H[i])^2)) / exp(β[i]) : 0.0 for i in 1:3]
    x1 = [i * gridinfo.h[1] for i in -Int(ceil(gridinfo.N_pad[1] / 2.0)):Int(ceil(gridinfo.N_pad[1] / 2.0)) - 1]
    x2 = [i * gridinfo.h[2] for i in -Int(ceil(gridinfo.N_pad[2] / 2.0)):Int(ceil(gridinfo.N_pad[2] / 2.0)) - 1]
    x3 = [i * gridinfo.h[3] for i in -Int(ceil(gridinfo.N_pad[3] / 2.0)):Int(ceil(gridinfo.N_pad[3] / 2.0)) - 1]

    DFT_fx = fft(fw[1].(x1))
    DFT_fy = fft(fw[2].(x2))
    DFT_fz = fft(fw[3].(x3))

    k1 = [2π * j / (gridinfo.L[1] + gridinfo.h[1] * 2 * gridinfo.pad[1]) for j in 0:(gridinfo.N_pad[1] - 1)]
    for j in 1:gridinfo.N_pad[1]
        if j - 1 > ceil(gridinfo.N_pad[1] / 2)
            k1[j] -= 2π * gridinfo.N_pad[1] / (gridinfo.L[1] + gridinfo.h[1] * 2 * gridinfo.pad[1])
        end
    end

    k2 = [2π * j / (gridinfo.L[2] + gridinfo.h[2] * 2 * gridinfo.pad[2]) for j in 0:(gridinfo.N_pad[2] - 1)]
    for j in 1:gridinfo.N_pad[2]
        if j - 1 > ceil(gridinfo.N_pad[2] / 2)
            k2[j] -= 2π * gridinfo.N_pad[2] / (gridinfo.L[2] + gridinfo.h[2] * 2 * gridinfo.pad[2])
        end
    end

    k3 = [2π * j / (gridinfo.L[3] + gridinfo.h[3] * 2 * gridinfo.pad[3]) for j in 0:(gridinfo.N_pad[3] - 1)]
    for j in 1:gridinfo.N_pad[3]
        if j - 1 > ceil(gridinfo.N_pad[3] / 2)
            k3[j] -= 2π * gridinfo.N_pad[3] / (gridinfo.L[3] + gridinfo.h[3] * 2 * gridinfo.pad[3])
        end
    end

    id_x = sortperm(k1)
    id_y = sortperm(k2)
    id_z = sortperm(k3)
    sort!(k1)
    sort!(k2)
    sort!(k3)

    DFT_fx_new = (gridinfo.L[1] + gridinfo.h[1] * 2 * gridinfo.pad[1]) / gridinfo.N_pad[1] .* DFT_fx[id_x]
    DFT_fy_new = (gridinfo.L[2] + gridinfo.h[2] * 2 * gridinfo.pad[2]) / gridinfo.N_pad[2] .* DFT_fy[id_y]
    DFT_fz_new = (gridinfo.L[3] + gridinfo.h[3] * 2 * gridinfo.pad[3]) / gridinfo.N_pad[3] .* DFT_fz[id_z]

    lif_x = linear_interpolation(k1,DFT_fx_new)
    lif_y = linear_interpolation(k2,DFT_fy_new)
    lif_z = linear_interpolation(k3,DFT_fz_new)
    F_fw = [x -> lif_x(x), x -> lif_y(x), x -> lif_z(x)]

    cheb_coefs = tuple([ChebCoef(fw[i], gridinfo.h[i], w[i], cheb_order) for i in 1:3]...)
    func_scale = (kx, ky, kz) -> (F_fw[1](kx) * F_fw[2](ky) * F_fw[3](kz))^(-2) * Tdk(kx, ky, kz, uspara, M_mid)
    scalefactor = ScalingFactor(func_scale, gridinfo)
    return (gridinfo, gridbox, cheb_coefs, scalefactor)
end


extra_pad_ratio = 8
cheb_order = 10
E_N_grid=[]
N_real_true = (40, 40, 40)
gridinfo_T, gridbox_T, cheb_coefs_T, scalefactor_T = mid_paras_gen_diffwindow(N_real_true, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid, alpha)
E_direct_true = energy_mid(qs, poses, gridinfo_T, gridbox_T, cheb_coefs_T, scalefactor_T)
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

        df = DataFrame(Window_type = 1, energy_exact = E_direct_true, abs_error = error, relative_error = error_rel, N_real = N_real)
        #CSV.write("data/Acc_T1_2_WinFunc_Gaussian.csv", df, append = true)
        CSV.write("data_for_manu/Acc_T3_2_Nspec_ES.csv", df, append = true)
end
