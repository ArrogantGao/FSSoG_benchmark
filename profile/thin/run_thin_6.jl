using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2, StatProfilerHTML, BenchmarkTools

Lxy0 = 100.0
Lz0 = 1.0

N_real0 = (160, 160)
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

for data in ["n_1000.jld2"]

    path = "../manuscript/reference/thin/$data" 
	@load path n_atoms Lx Lz atoms info energy_ewald
	boundary = ExTinyMD.Q2dBoundary(Lx, Lx, Lz)

    ratio = (n_atoms / 1000)^(1/2)
	N_real = Int.(ceil.(ratio .* N_real0))

	fssog_interaction = FSSoGThinInteraction((Lx, Lx, Lz), n_atoms, r_c, Q, 0.5, N_real, R_z, w, β, cheb_order, Taylor_Q, Rz_0, Q_0; preset = preset, ϵ = 1.0)
	fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
    # warm up
    runtime = @belapsed ExTinyMD.energy($fssog_interaction, $fssog_neighbor, $info, $atoms)
    @show runtime
    
    # run
    @info "profiling $n_atoms"
    @profilehtml for i in 1:100 ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms) end
end
