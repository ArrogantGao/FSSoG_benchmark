using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

L0 = 1.0
N_real0 = (100, 100, 100)

w = (22, 22, 22)
β = 5.0 .* w
target = 1.8
extra_pad_ratio_intial = ((target-1.0)*N_real0[1]+6.8*(target-1.0)*w[1]+2.4*target-4.4)/(2*w[1]) 
@show extra_pad_ratio_intial
cheb_order = 10
preset = 6
uspara = USeriesPara(preset)
M_mid_initial = 12
eta = uspara.sw[M_mid_initial][1] / L0 + 0.0001
@show eta
N_grid = (3, 3, 18)
Qs = 32
Ql = 20
r_c = 0.499
Q_0 = 32
Rz_0 = 25


for data in ["n_100.jld2","n_318.jld2","n_1000.jld2", "n_3164.jld2", "n_10000.jld2"]

    path = "./reference/cube/$data" 
	@load path n_atoms L atoms info energy_ewald 

	@show n_atoms, L, energy_ewald

	boundary = ExTinyMD.Q2dBoundary(L, L, L)

    ratio = (n_atoms / 100)^(1/3)
	N_real = Int.(ceil.(ratio .* N_real0))
	extra_pad_ratio = Int(ceil((extra_pad_ratio_intial.*ratio)))
	#extra_pad_ratio = Int(ceil(((target-1.0)*N_real[1]+6.8*(target-1.0)*w[1]+2.4*target-4.4)/(2*w[1])))

	@show N_real, extra_pad_ratio

	M_mid = proper_M(eta, L, uspara)
	@show M_mid

	@info "FSSoG, n_atoms = $n_atoms"
	fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Qs, 0.1, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Ql, Rz_0, Q_0; preset = preset, ϵ = 1.0)
	fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
	@time energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

	λ = fssog_interaction.gridinfo.N_pad[3] / fssog_interaction.gridinfo.N_image[3]
	@show λ


	abs_error = abs(energy_ewald - energy_sog)
	relative_error = abs(energy_ewald - energy_sog) / abs(energy_ewald)

	@show energy_sog, abs_error, relative_error

	# abs_error = abs(energy_ewald - energy_sog)
	# relative_error = abs(energy_ewald - energy_sog) / abs(energy_ewald)

	# @show energy_sog, abs_error, relative_error

	df = DataFrame(n_atoms = n_atoms, E_exact = energy_ewald, E_fssog = energy_sog, abs_error = abs_error, relative_error = relative_error)
	CSV.write("data/Acc_T6_cube_total_6.csv", df, append = true)
end
