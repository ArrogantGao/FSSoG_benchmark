using FastSpecSoG, BenchmarkTools, ExTinyMD

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
position = [info.particle_info[i].position for i in 1:n_atoms]
charge = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]

# here 10.0 is r_c, 3.0 is k_c (not important here)
FastSpecSoG_interaction = FSSoG_naive((L, L, L), n_atoms, 10.0, 3.0, preset = 5)
FastSpecSoG_neighbor = CellList3D(info, FastSpecSoG_interaction.r_c, boundary, 1)
Es_naive = short_energy_naive(FastSpecSoG_interaction, FastSpecSoG_neighbor, position, charge)

@benchmark short_energy_naive($FastSpecSoG_interaction, $FastSpecSoG_neighbor, $position, $charge)