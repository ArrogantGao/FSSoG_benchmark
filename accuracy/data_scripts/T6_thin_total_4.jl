using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

function main()

	N_real0 = (24, 24)
	R_z = 4
	w = (6, 6)
	β = 7.5 .* w
	cheb_order = 8
	preset = 4
	Q = 16
	Q_0 = 8
	Rz_0 = 6
	Taylor_Q = 4
	r_c = 10.0

	data_file = "data/Acc_T6_thin_total_4.csv"
	CSV.write(data_file, DataFrame(n_atoms = Int[], abs_error = Float64[], relative_error = Float64[]))

	for data in ["n_1000.jld2", "n_3164.jld2", "n_10000.jld2", "n_31624.jld2", "n_100000.jld2", "n_316228.jld2", "n_1000000.jld2"]

		path = "./reference/thin_200/$data" 
		@load path n_atoms Lx Lz atoms info energy_per_atoms

		@show n_atoms, Lx, Lz

		boundary = ExTinyMD.Q2dBoundary(Lx, Lx, Lz)

		ratio = (n_atoms / 1000)^(1/2)
		N_real = Int.(ceil.(ratio .* N_real0))

		@info "FSSoG, n_atoms = $n_atoms"
		fssog_interaction = FSSoGThinInteraction((Lx, Lx, Lz), n_atoms, r_c, Q, 0.001, N_real, R_z, w, β, cheb_order, Taylor_Q, Rz_0, Q_0; preset = preset, ϵ = 1.0)
		fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
		@time energy_sog_per_atom = energy_per_atom(fssog_interaction, fssog_neighbor, info, atoms)

		abs_error = rmsd(energy_sog_per_atom[1:length(energy_per_atoms)], energy_per_atoms)
		relative_error = rrmsd(energy_sog_per_atom[1:length(energy_per_atoms)], energy_per_atoms)

		@show abs_error, relative_error

		df = DataFrame(n_atoms = n_atoms, abs_error = abs_error, relative_error = relative_error)
		CSV.write(data_file, df, append = true)
	end
end

main()