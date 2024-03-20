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
energy_short_Q = [[] for i in 1:5]
cheb_orders = [1:3:40...]

# fix r_c = 10.0
r_c = 10.0

for preset in 1:5
    exact_interaction = FSSoG_naive((L, L, L), n_atoms, 49.9, 3.0, preset = preset)
    exact_neighbor = CellList3D(info, exact_interaction.r_c, boundary, 1)
    Es_exact = short_energy_naive(exact_interaction, exact_neighbor, position, charge)
    @show preset, Es_exact
    for Q in cheb_orders
        uspara_cheb, F0 = Es_Cheb_precompute(preset, 0.5, 10.0, Q)
        neighbor = CellList3D(info, r_c, boundary, 1)

        Es_cheb = short_energy_Cheb(uspara_cheb, 10.0, F0, boundary, neighbor, position, charge)
        push!(energy_short_Q[preset], Es_cheb)
    end
    push!(energy_short_exact, Es_exact)
end

fig_abs = plot(dpi = 500, xlabel = "Cheb Order", ylabel = "absolute error", title = L"r_c = 10", ylim = [1e-12, 10.0])
for i in 1:5
    plot!(fig_abs, cheb_orders, abs.(energy_short_Q[i] .- energy_short_exact[i]), label = "preset = $i", yscale = :log10, marker = :circle)
end

savefig(fig_abs, "figs/accuracy_short_cheb_total.png")