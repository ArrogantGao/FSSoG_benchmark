using FastSpecSoG, ExTinyMD, Plots, LaTeXStrings, CSV, DataFrames

n_array = [8, 16, 32, 64, 128]

errors = [[] for i in 1:2]

for preset in 1:2
    M_mid = 0
    n_atoms = 100
    L = (100.0, 100.0, 1.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    uspara = USeriesPara(preset)
    # k = π N / L, exp(- s_min^2 k^2 / 4) = 1e-16
    sm = uspara.sw[1][1]
    N = 50
    E_direct_k = long_energy_us_k(qs, poses, N, L, uspara, M_mid + 1, length(uspara.sw))

    for n in n_array
        N_real = (n, n)
        R_z = 16
        w = N_real .÷ 8
        β = 5.0 .* w
        cheb_order = 16
        Taylor_Q = 16

        gridinfo, pad_grids, cheb_coefs, scalefactors, H_r, H_c, cheb_value, r_z = thin_paras_gen(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)

        E_thin_k = energy_long_thin_k(qs, poses, L, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value)

        df = DataFrame(preset = preset, Nxy = n, E_thin_k = E_thin_k, E_direct_k = E_direct_k, error = abs(E_thin_k - E_direct_k) / abs(E_direct_k))
        CSV.write("data/thin_long_Nxy.csv", df, append = true)
        push!(errors[preset], abs(E_thin_k - E_direct_k) / abs(E_direct_k))
        @show preset, n, abs(E_thin_k - E_direct_k) / abs(E_direct_k)
    end
end

df = CSV.read("data/thin_long_Nxy.csv", DataFrame)
fig = plot(xlabel = "log2(Nxy)", ylabel = "Relative error", title = "Long range energy (thin)", dpi = 500)
for i in 1:2
    dfi = filter(row -> row.preset == i, df)
    plot!(log2.(dfi.Nxy), log10.(dfi.error), label = "set: $i", ylims = (-16, 1), marker = :circle)
end
savefig(fig, "figs/long_range_thin_errors_Nxy.png")