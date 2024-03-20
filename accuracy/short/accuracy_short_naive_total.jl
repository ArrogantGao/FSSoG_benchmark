using FastSpecSoG, BenchmarkTools, ExTinyMD, Plots, LaTeXStrings
using Random
Random.seed!(1234)

n_atoms = 1000
L = 100.0
boundary = ExTinyMD.Q2dBoundary(L, L, L)

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms÷2
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
end

for i in n_atoms÷2 + 1 : n_atoms
    push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)
position = [info.particle_info[i].position.coo for i in 1:n_atoms]
charge = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]

energy_short_exact = []
energy_short_rc = [[] for i in 1:5]
r_c = [1.0:0.5:15.0...]
for preset in 1:5
    exact_interaction = FSSoG_naive((L, L, L), n_atoms, 49.9, 3.0, preset = preset)
    exact_neighbor = CellList3D(info, exact_interaction.r_c, boundary, 1)
    Es_exact = short_energy_naive(exact_interaction, exact_neighbor, position, charge)
    @show preset, Es_exact
    for rc in r_c
        FastSpecSoG_interaction = FSSoG_naive((L, L, L), n_atoms, rc, 3.0, preset = preset)
        FastSpecSoG_neighbor = CellList3D(info, FastSpecSoG_interaction.r_c, boundary, 1)

        Es_naive = short_energy_naive(FastSpecSoG_interaction, FastSpecSoG_neighbor, position, charge)
        push!(energy_short_rc[preset], Es_naive)
    end
    push!(energy_short_exact, Es_exact)
end

fig = plot(dpi = 500, xlabel = "r_c", ylabel = "relative error", ylim = [1e-12, 10.0])
for i in 1:5
    plot!(fig, r_c, abs.(energy_short_rc[i] .- energy_short_exact[i]) ./ abs(energy_short_exact[i]), label = "preset = $i", yscale = :log10, marker = :circle)
end

savefig(fig, "figs/accuracy_short_naive_relative.png")

fig_abs = plot(dpi = 500, xlabel = L"r_c", ylabel = "absolute error", ylim = [1e-12, 10.0])
for i in 1:5
    plot!(fig_abs, r_c, abs.(energy_short_rc[i] .- energy_short_exact[i]), label = "preset = $i", yscale = :log10, marker = :circle)
end

savefig(fig_abs, "figs/accuracy_short_naive_absolute.png")