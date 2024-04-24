using EwaldSummations, ExTinyMD, JLD2, FastSpecSoG
using Random
Random.seed!(123)

n_atoms0 = 100
Lx0 = 1.0
Lz = 0.1

s = 6.0
ratios = [1000.0]

for ratio in ratios
    n_atoms = Int(ceil(n_atoms0 * ratio))
    if isodd(n_atoms)
        n_atoms += 1
    end

    Lx = Lx0 * ratio^(0.5)
    
    boundary = ExTinyMD.Q2dBoundary(Lx, Lx, Lz)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, Lx, 0.0, Lx, 0.0, Lz), boundary; min_r = 0.05, temp = 1.0)

    # s / α = Lx / 2 - 0.1
    α = s / (Lx / 2 - 0.1)
    Ewald2D_interaction = Ewald2DInteraction(n_atoms, s, α, (Lx, Lx, Lz), ϵ = 1.0)
    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)

    r_c = Ewald2D_interaction.r_c

    @show ratio, n_atoms, Lx, Lz, n_atoms/(Lx^2 * Lz), α, r_c

    charge = [atoms[info.particle_info[i].id].charge for i in 1:Ewald2D_interaction.n_atoms]
    position = [info.particle_info[i].position for i in 1:Ewald2D_interaction.n_atoms]

    # @time energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)
    @time energy_long = Ewald2D_long_energy(Ewald2D_interaction, position, charge)
    @time energy_short = Ewald2D_short_energy(Ewald2D_interaction, Ewald2D_neighbor, position, charge)
    energy_ewald = energy_long + energy_short
    @show energy_ewald


	# N_real = Int.(ceil.(ratio^(1/2) .* N_real0))

	# @info "FSSoG, n_atoms = $n_atoms"
	# fssog_interaction = FSSoGThinInteraction((Lx, Lx, Lz), n_atoms, r_c, Q, 0.3, N_real, R_z, w, β, cheb_order, Taylor_Q, Rz_0, Q_0; preset = preset, ϵ = 1.0)
	# fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
	# energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)
    # @show energy_sog

    @save "reference/thin/high_rho/n_$(n_atoms).jld2" n_atoms Lx Lz atoms info energy_ewald
end
