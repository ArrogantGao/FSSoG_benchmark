using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

Lxy0 = 100.0
Lz0 = 1.0

N_real0 = (64, 64)
R_z = 6
w = (8, 8)
β = 7.5 .* w
cheb_order = 4
preset = 2
Q = 4
Q_0 = 4
Taylor_Q = 4
r_c = 10.0
Rz_0 = 4

for data in ["n_1000.jld2", "n_3164.jld2", "n_10000.jld2", "n_31624.jld2", "n_100000.jld2"]

    path = "./reference/thin/$data" 
	@load path n_atoms Lx Lz atoms info energy_ewald 

	@show n_atoms, Lx, Lz, energy_ewald

	boundary = ExTinyMD.Q2dBoundary(Lx, Lx, Lz)

    ratio = (n_atoms / 1000)^(1/2)
	N_real = Int.(ceil.(ratio .* N_real0))

	@info "FSSoG, n_atoms = $n_atoms"
	fssog_interaction = FSSoGThinInteraction((Lx, Lx, Lz), n_atoms, r_c, Q, 0.5, N_real, R_z, w, β, cheb_order, Taylor_Q, Rz_0, Q_0; preset = preset, ϵ = 1.0)
	fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
	@time energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

	abs_error = abs(energy_ewald - energy_sog)
	relative_error = abs(energy_ewald - energy_sog) / abs(energy_ewald)

	@show energy_sog, abs_error, relative_error

	df = DataFrame(n_atoms = n_atoms, E_exact = energy_ewald, E_fssog = energy_sog, abs_error = abs_error, relative_error = relative_error)
	CSV.write("data/Acc_T6_thin_total_2.csv", df, append = true)
end
