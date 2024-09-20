using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

function main()

	L0 = 20.0

	N_real0 = (16, 16, 16)
	w = (6, 6, 6)
	β = 5.0 .* w
	extra_pad_ratio_intial = 3
	lambda_true = (2*extra_pad_ratio_intial*w[1]+N_real0[1]+6.8*w[1]+4.4)/(N_real0[1]+6.8*w[1]+2.4)
	@show lambda_true
	cheb_order = 8
	preset = 4
	uspara = USeriesPara(preset) 
	M_mid_initial = 6
	eta = uspara.sw[M_mid_initial][1] / L0 + 0.0001
	@show eta
	N_grid = (1, 1, 9)
	Qs = 16
	Ql = 8
	r_c = 9.99
	Rz_0 = 9
	Q_0 = 8

	CSV.write("data/Acc_T6_cube_total_4.csv", DataFrame(n_atoms = [], abs_error = [], relative_error = []))

	for data in ["n_1000.jld2", "n_3164.jld2", "n_10000.jld2", "n_31624.jld2", "n_100000.jld2", "n_316228.jld2", "n_1000000.jld2"]

		path = "./reference/cube_200/$data" 
		@load path n_atoms L atoms info energy_per_atoms

		@show n_atoms, L

		boundary = ExTinyMD.Q2dBoundary(L, L, L)

		ratio = (n_atoms / 1000)^(1/3)
		N_real = Int.(ceil.(ratio .* N_real0))
		extra_pad_ratio = Int(ceil(extra_pad_ratio_intial * ratio))

		M_mid = proper_M(eta, L, uspara)
		@show M_mid

		@info "FSSoG, n_atoms = $n_atoms"
		fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Qs, 0.01, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Ql, Rz_0, Q_0; preset = preset, ϵ = 1.0)
		fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
		@time energy_sog_per_atom = energy_per_atom(fssog_interaction, fssog_neighbor, info, atoms)

		abs_error = rmsd(energy_sog_per_atom[1:length(energy_per_atoms)], energy_per_atoms)
		relative_error = rrmsd(energy_sog_per_atom[1:length(energy_per_atoms)], energy_per_atoms)

		@show abs_error, relative_error
		
		df = DataFrame(n_atoms = n_atoms, abs_error = abs_error, relative_error = relative_error)
		CSV.write("data/Acc_T6_cube_total_4.csv", df, append = true)
	end
end

main()