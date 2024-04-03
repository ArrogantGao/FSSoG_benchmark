using FastSpecSoG, ExTinyMD, Plots, LaTeXStrings, CSV, DataFrames

<<<<<<< HEAD
n_array = [1:1:8...]

pm = 3
Ns = [30, 50, 100]
=======
n_array = [1:1:15...]

pm = 6
>>>>>>> main
errors = [[] for i in 1:pm]

for preset in 1:pm
    n_atoms = 100
    L = (100.0, 100.0, 1.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    uspara = USeriesPara(preset)
<<<<<<< HEAD
    E_direct_k = long_energy_us_k(qs, poses, Ns[preset], L, uspara, 1, length(uspara.sw))

    for n in n_array
        N_real = (128, 128)
        R_z = 16
        w = N_real .÷ 4
        β = 5.0 .* w
=======
    E_direct_k = long_energy_us_k(qs, poses, 1e-18, L, uspara, 1, length(uspara.sw))

    for n in n_array
        N_real = (256, 256)
        R_z = 16
        w = N_real .÷ 4
        β = 8.0 .* w
>>>>>>> main
        cheb_order = 16
        Taylor_Q = n

        gridinfo, pad_grids, cheb_coefs, scalefactors, H_r, H_c, cheb_value, r_z = thin_paras_gen(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)

        E_thin_k = energy_long_thin_k(qs, poses, L, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value)

        df = DataFrame(preset = preset, Taylor_Q = n, E_thin_k = E_thin_k, E_direct_k = E_direct_k, error = abs(E_thin_k - E_direct_k) / abs(E_direct_k))
        CSV.write("data/thin_long_taylorQ.csv", df, append = true)
        push!(errors[preset], abs(E_thin_k - E_direct_k) / abs(E_direct_k))
        @show preset, n, abs(E_thin_k - E_direct_k) / abs(E_direct_k)
    end
end

fig = plot(xlabel = "Taylor Q", ylabel = "Relative error", title = "Long range energy (thin)", dpi = 500)
for i in 1:pm
    plot!(n_array, log10.(errors[i]), label = "set: $i", ylims = (-16, 1), marker = :circle)
end
savefig(fig, "figs/long_range_thin_errors_taylorQ.png")