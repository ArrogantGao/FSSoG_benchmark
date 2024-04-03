using FastSpecSoG, ExTinyMD, Plots, LaTeXStrings, CSV, DataFrames

n_array = [8, 16, 32, 64, 128, 256]

pm = 6
errors = [[] for i in 1:pm]
M_mid = 0

for preset in 1:pm
    n_atoms = 100
    L = (100.0, 100.0, 1.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    uspara = USeriesPara(preset)
    E_direct_k = long_energy_us_k(qs, poses, 1e-18, L, uspara, 1, length(uspara.sw))

    for n in n_array

        N_grid = (n, n, 16)

        k_x, k_y, r_z, H_r, H_c, phase_x, phase_y, sort_z, z = long_paras_gen(L, N_grid, n_atoms)
        cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, 16)
        E_thin_k = energy_long_cheb_k(qs, poses, L, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, cheb_mat)

        df = DataFrame(preset = preset, Nxy = n, E_thin_k = E_thin_k, E_direct_k = E_direct_k, error = abs(E_thin_k - E_direct_k) / abs(E_direct_k))
        CSV.write("data/thin_long_Nxy_cheb.csv", df, append = true)
        push!(errors[preset], abs(E_thin_k - E_direct_k) / abs(E_direct_k))
        @show preset, n, abs(E_thin_k - E_direct_k) / abs(E_direct_k)
    end
end

fig = plot(xlabel = "Nxy", ylabel = "Relative error", title = "Long range energy (thin cheb)", dpi = 500)
for i in 1:pm
    plot!(log2.(n_array), log10.(errors[i]), label = "set: $i", ylims = (-16, 1), marker = :circle)
end
savefig(fig, "figs/long_range_thin_errors_Nxy_cheb.png")