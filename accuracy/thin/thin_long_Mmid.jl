using FastSpecSoG, ExTinyMD, Plots, LaTeXStrings, CSV, DataFrames

n_array = [8, 16, 32, 64]

errors = [[] for i in 1:6]
Mmid_array = [0, 1, 2, 3, 4, 5]
preset = 2

for M_mid in Mmid_array
    n_atoms = 100
    L = (100.0, 100.0, 1.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    uspara = USeriesPara(preset)
    E_direct_k = long_energy_us_k(qs, poses, 50, L, uspara, M_mid + 1, length(uspara.sw))

    for n in n_array

        N_grid = (n, n, 16)
        k_x, k_y, r_z, H_r, H_c, phase_x, phase_y,  sort_z, z = long_paras_gen(L, N_grid, n_atoms)
        cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, 64)
        E_thin_k = energy_long_cheb_k(qs, poses, L, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, cheb_mat)

        df = DataFrame(M_mid = M_mid, Nxy = n, E_thin_k = E_thin_k, E_direct_k = E_direct_k, error = abs(E_thin_k - E_direct_k) / abs(E_direct_k))
        CSV.write("data/thin_long_Mmid.csv", df, append = true)
        push!(errors[preset], abs(E_thin_k - E_direct_k) / abs(E_direct_k))
        @show M_mid, n, abs(E_thin_k - E_direct_k) / abs(E_direct_k)
    end
end

df = CSV.read("data/thin_long_Mmid.csv", DataFrame)
fig = plot(xlabel = "log2(Nxy)", ylabel = "Relative error", title = "Long range energy (thin)", dpi = 500)
for i in 0:5
    dfi = filter(row -> row.M_mid == i, df)
    plot!(log2.(dfi.Nxy), log10.(dfi.error), label = "M_mid = $i", ylims = (-16, 1), marker = :circle)
end
savefig(fig, "figs/long_range_thin_errors_Mmid.png")