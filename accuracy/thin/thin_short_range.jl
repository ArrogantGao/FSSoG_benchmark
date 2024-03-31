using FastSpecSoG, CSV, DataFrames, Plots, ExTinyMD, LaTeXStrings
using Random
Random.seed!(1234)


n_atoms = 100
Lxy = 100.0
Lz = 1.0
boundary = ExTinyMD.Q2dBoundary(Lxy, Lxy, Lz)

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms÷2
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
end

for i in n_atoms÷2 + 1 : n_atoms
    push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
end

r_cs = [0.5:0.5:20.0...]

info = SimulationInfo(n_atoms, atoms, (0.0, Lxy, 0.0, Lxy, 0.0, Lz), boundary; min_r = 1.0, temp = 1.0)
position = [info.particle_info[i].position.coo for i in 1:n_atoms]
charge = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]

energy_short_exact = []
energy_short_rc = [[] for i in 1:6]

Q = 64

for preset in 1:6
    exact_interaction = FSSoG_naive((Lxy, Lxy, Lz), n_atoms, 49.9, 3.0, preset = preset)
    exact_neighbor = CellList3D(info, exact_interaction.r_c, boundary, 1)
    Es_exact = short_energy_naive(exact_interaction, exact_neighbor, position, charge)
    push!(energy_short_exact, Es_exact)
    @show preset, Es_exact

    for r_c in r_cs
        uspara_cheb, F0 = Es_Cheb_precompute(preset, 0.5, r_c, Q)
        neighbor = CellList3D(info, r_c, boundary, 1)
        Es_cheb = short_energy_Cheb(uspara_cheb, r_c, F0, boundary, neighbor, position, charge)
        push!(energy_short_rc[preset], Es_cheb)
    end
end

errors = [abs.((energy_short_rc[i] .- energy_short_exact[i]) ./ energy_short_exact[i]) for i in 1:6]

fig = plot(xlabel = "r_c", ylabel = "Relative error", title = "Short range energy (thin)", dpi = 500)
for i in 1:6
    plot!(r_cs, log10.(errors[i]), label = "set: $i", ylims = (-16, 1), marker = :circle)
end
savefig(fig, "figs/short_range_thin_errors.png")

# df = DataFrame(r_c = r_cs, errors...)
# CSV.write("data/short_range_thin_errors.csv", df)