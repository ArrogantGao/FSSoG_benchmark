using EwaldSummations, ExTinyMD, JLD2, FastSpecSoG
using Random
Random.seed!(123)

n_atoms0 = 1000
Lx0 = 30.0
Lz = 0.3

s = 7.0
ratios = [10^i for i in 0.0:0.5:2.0]

# N_real0 = (160, 160)
N_real0 = (96, 96)
R_z = 12
w = (16, 16)
β = 7.5 .* w
cheb_order = 12
preset = 6
Q = 32
Q_0 = 16
Rz_0 = 12
Taylor_Q = 8
r_c = 10.0

for ratio in ratios
    n_atoms = Int(ceil(n_atoms0 * ratio))
    if isodd(n_atoms)
        n_atoms += 1
    end

    Lx = Lx0 * ratio^(1/2)
    
    boundary = ExTinyMD.Q2dBoundary(Lx, Lx, Lz)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, Lx, 0.0, Lx, 0.0, Lz), boundary; min_r = 0.3, temp = 1.0)

    #s / α = Lx / 2 - 0.1
    # α = s / (Lx0 / 2 - 0.1) * Lx0 / Lx
    # Ewald2D_interaction = Ewald2DInteraction(n_atoms, s, α, (Lx, Lx, Lz), ϵ = 1.0)

    # r_c = Ewald2D_interaction.r_c

    # @show ratio, n_atoms, Lx, Lz, n_atoms/(Lx^2 * Lz), α, r_c

    # Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    # @time energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)
    # @show energy_ewald
    energy_ewald = 0.0


	N_real = Int.(ceil.(ratio .* N_real0))

	@info "FSSoG, n_atoms = $n_atoms"
	fssog_interaction = FSSoGThinInteraction((Lx, Lx, Lz), n_atoms, r_c, Q, 0.3, N_real, R_z, w, β, cheb_order, Taylor_Q, Rz_0, Q_0; preset = preset, ϵ = 1.0)
	fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
	energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)
    @show energy_sog

    @save "reference/thin/n_$(n_atoms).jld2" n_atoms Lx Lz atoms info energy_ewald energy_sog
end
