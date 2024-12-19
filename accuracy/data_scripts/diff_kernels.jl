function Tdk(kx::T, ky::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    k2 = kx^2 + ky^2

    Tdk_value = zero(T)
    for i in 1:M_mid
        sl, wl = uspara.sw[i]
        Tdk_value += wl * sl^3 * exp(- sl^2 * k2 / 4)
    end

    return sqrt(π)/4 * Tdk_value
end

function thin_grids_gen_ES(N_real::NTuple{2, Int}, R_z::Int, w::NTuple{2, Int}, β::NTuple{2, T}, L::NTuple{3, T}, cheb_order::Int, uspara::USeriesPara{T}, Taylor_Q::Int) where{T}
    periodicity = (true, true)
    extra_pad = (0, 0)

    gridinfo = GridInfo(N_real, w, periodicity, extra_pad, (L[1], L[2]))
    pad_grids = [zeros(Complex{T}, N_real[1], N_real[2], R_z) for i in 1:Taylor_Q]

    k_x = gridinfo.k[1]
    k_y = gridinfo.k[2]
    taylor_mats = TaylorUSeries(k_x, k_y, L[3], uspara, 0, Taylor_Q)

    H = [(w[i] + 0.5) * gridinfo.h[i] for i in 1:2]
    fw = [x -> abs(x) < H[i] ? exp(β[i] * sqrt(1 - (x / H[i])^2)) / exp(β[i]) : 0.0 for i in 1:2]
    cheb_coefs = tuple([ChebCoef(fw[i], gridinfo.h[i], w[i], cheb_order) for i in 1:2]...)

    x1 = [i * gridinfo.h[1] for i in -Int(ceil(gridinfo.N_pad[1] / 2.0)):Int(ceil(gridinfo.N_pad[1] / 2.0)) - 1]
    x2 = [i * gridinfo.h[2] for i in -Int(ceil(gridinfo.N_pad[2] / 2.0)):Int(ceil(gridinfo.N_pad[2] / 2.0)) - 1]

    DFT_fx = fft(fw[1].(x1))
    DFT_fy = fft(fw[2].(x2))

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

    id_x = sortperm(k1)
    id_y = sortperm(k2)
    sort!(k1)
    sort!(k2)

    DFT_fx_new = (gridinfo.L[1] + gridinfo.h[1] * 2 * gridinfo.pad[1]) / gridinfo.N_pad[1] .* DFT_fx[id_x]
    DFT_fy_new = (gridinfo.L[2] + gridinfo.h[2] * 2 * gridinfo.pad[2]) / gridinfo.N_pad[2] .* DFT_fy[id_y]

    lif_x = linear_interpolation(k1,DFT_fx_new)
    lif_y = linear_interpolation(k2,DFT_fy_new)
    F_fw = [x -> lif_x(x), x -> lif_y(x)]

    func_scale = (k_x, k_y) -> (F_fw[1](k_x) * F_fw[2](k_y))^(-2)
    sf0 = ScalingFactor(func_scale, gridinfo)
    scalefactors = [ScalingFactor(func_scale, sf0.factors .* taylor_mats[i]) for i in 1:Taylor_Q]

    return (gridinfo, pad_grids, cheb_coefs, scalefactors)
end

function thin_grids_gen_Gaussian(N_real::NTuple{2, Int}, R_z::Int, w::NTuple{2, Int}, alpha::NTuple{2, T}, L::NTuple{3, T}, cheb_order::Int, uspara::USeriesPara{T}, Taylor_Q::Int) where{T}
    periodicity = (true, true)
    extra_pad = (0, 0)

    gridinfo = GridInfo(N_real, w, periodicity, extra_pad, (L[1], L[2]))
    pad_grids = [zeros(Complex{T}, N_real[1], N_real[2], R_z) for i in 1:Taylor_Q]

    k_x = gridinfo.k[1]
    k_y = gridinfo.k[2]
    taylor_mats = TaylorUSeries(k_x, k_y, L[3], uspara, 0, Taylor_Q)
    
    H = [(w[i] + 0.5) * gridinfo.h[i] for i in 1:2]
    f_window = [x -> abs(x) < H[i] ? exp(-alpha[i] * (x/H[i])^2) : 0.0 for i in 1:2]
    F_f_window = [x -> sqrt(π/alpha[i]) * H[i] * exp(-x^2 * H[i]^2 / (4 * alpha[i]) )  for i in 1:2]
    cheb_coefs = tuple([ChebCoef(f_window[i], gridinfo.h[i], w[i], cheb_order) for i in 1:2]...)

    k_x = gridinfo.k[1]
    k_y = gridinfo.k[2]
    taylor_mats = TaylorUSeries(k_x, k_y, L[3], uspara, 0, Taylor_Q)

    func_scale = (k_x, k_y) -> (F_f_window[1](k_x) * F_f_window[2](k_y))^(-2)
    sf0 = ScalingFactor(func_scale, gridinfo)
    scalefactors = [ScalingFactor(func_scale, sf0.factors .* taylor_mats[i]) for i in 1:Taylor_Q]

    return (gridinfo, pad_grids, cheb_coefs, scalefactors)
end
