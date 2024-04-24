using FastSpecSoG, ExTinyMD, BenchmarkTools, Plots, CSV, DataFrames, ArgParse, JLD2

function parse_commandline()

	s = ArgParseSettings()

	@add_arg_table s begin
		"--output-file"
		help = "relative file path of the output file"
	end

	return parse_args(s)

end

function main()

	parsed_args = parse_commandline()
	output_file = parsed_args["output-file"]

	N_real0 = (320, 320)
	R_z = 8
	w = (7, 7)
	β = 7.5 .* w
	cheb_order = 8
	preset = 6
	Q = 32
	Q_0 = 8
	Rz_0 = 8
	Taylor_Q = 6
	r_c = 1.799918355128706

	for data in ["n_100000.jld2"]

		path = "../manuscript/reference/thin/high_rho/$data" 
		@load path n_atoms Lx Lz atoms info energy_ewald

		boundary = ExTinyMD.Q2dBoundary(Lx, Lx, Lz)

		ratio = 1.0
		N_real = Int.(ceil.(ratio .* N_real0))
		interaction = FSSoGThinInteraction((Lx, Lx, Lz), n_atoms, r_c, Q, 0.05, N_real, R_z, w, β, cheb_order, Taylor_Q, Rz_0, Q_0; preset = preset, ϵ = 1.0)
		neighbor = CellList3D(info, interaction.r_c, boundary, 1)

		for i in 1:interaction.n_atoms
			interaction.charge[i] = atoms[info.particle_info[i].id].charge
			interaction.position[i] = info.particle_info[i].position.coo
		end

		t_short = @belapsed energy_short($interaction, $neighbor)
		t_long = @belapsed energy_long($interaction)

		@show preset, n_atoms, t_short, t_long

		df = DataFrame(preset = preset, n_atoms = n_atoms, t_short = t_short, t_long = t_long, t_total = t_short + t_long)
		CSV.write(output_file, df, append = true)
	end
end

main()
