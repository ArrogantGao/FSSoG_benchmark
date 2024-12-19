using EwaldSummations, ExTinyMD, JLD2
using Random
Random.seed!(123)

n_atoms0 = 100
L0 = 1.0

s = 7.0
ratios = [1.0]

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

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 0.01, temp = 1.0)

    #s / α = L / 2 - 0.1
    α = s / (L / 2 - 0.1)
    interaction = Ewald2DInteraction(n_atoms, s, α, (L, L, L), ϵ = 1.0)

    r_c = interaction.r_c

    @show ratio, n_atoms, L, n_atoms/L^3, α, r_c

    neighbor = CellList3D(info, interaction.r_c, boundary, 1)

    fs = force(interaction, neighbor, info, atoms)
    force_ewald = [f.coo for f in fs]

    @save "reference/cube_L1/n_$(n_atoms)_force.jld2" n_atoms L atoms info force_ewald
end
