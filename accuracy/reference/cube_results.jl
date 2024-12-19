using EwaldSummations, ExTinyMD, JLD2
using Random
Random.seed!(1234)

n_atoms0 = 100
L0 = 1.0

s = 6.0
ratios = [10^i for i in 0.0:0.5:3.0]

for ratio in ratios
    n_atoms = Int(ceil(n_atoms0 * ratio))
    if isodd(n_atoms)
        n_atoms += 1
    end

    L = L0 * ratio^(1/3)
    
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 0.1, temp = 1.0)

    α = s / (L / 2 - 0.01)
    Ewald2D_interaction = Ewald2DInteraction(n_atoms, s, α, (L, L, L), ϵ = 1.0)

    r_c = Ewald2D_interaction.r_c
    k_c = Ewald2D_interaction.k_c

    @show ratio, n_atoms, L, n_atoms/L^3, α, r_c, k_c

    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)

    charge = [atoms[info.particle_info[i].id].charge for i in 1:Ewald2D_interaction.n_atoms]
    position = [info.particle_info[i].position for i in 1:Ewald2D_interaction.n_atoms]

    # @time energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)
    @time energy_long = Ewald2D_long_energy(Ewald2D_interaction, position, charge)
    @time energy_short = Ewald2D_short_energy(Ewald2D_interaction, Ewald2D_neighbor, position, charge)
    energy_ewald = energy_long + energy_short
    @show energy_ewald

    @save "reference/cube/n_$(n_atoms).jld2" n_atoms L atoms info energy_ewald
end
