using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

L0 = 20.0
N_real0 = (36, 36, 36)

w = (9, 9, 9)
β = 5.0 .* w
extra_pad_ratio_intial = 5.095
cheb_order = 10
preset = 6
uspara = USeriesPara(preset)
M_mid_initial = 14
eta = uspara.sw[M_mid_initial][1] / L0 + 0.0001
@show eta
N_grid = (2, 2, 15)
Qs = 28
Ql = 16
r_c = 9.99
Q_0 = 32
Rz_0 = 25

target = 2.0

for data in ["n_1000.jld2", "n_3164.jld2", "n_10000.jld2", "n_31624.jld2", "n_100000.jld2"]

    path = "./reference/cube/$data" 
	@load path n_atoms L atoms info energy_ewald energy_fssog

	@show n_atoms, L, energy_ewald, energy_fssog

	boundary = ExTinyMD.Q2dBoundary(L, L, L)

    ratio = (n_atoms / 1000)^(1/3)
	N_real = Int.(ceil.(ratio .* N_real0))
	extra_pad_ratio = Int(ceil((extra_pad_ratio_intial.*ratio)))

	@show N_real, extra_pad_ratio

	M_mid = proper_M(eta, L, uspara)
	@show M_mid

	@info "FSSoG, n_atoms = $n_atoms"
	fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Qs, 0.1, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Ql, Rz_0, Q_0; preset = preset, ϵ = 1.0)
	fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
	@time energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

	λ = fssog_interaction.gridinfo.N_pad[3] / fssog_interaction.gridinfo.N_image[3]
	@show λ


	abs_error = abs(energy_fssog - energy_sog)
	relative_error = abs(energy_fssog - energy_sog) / abs(energy_fssog)

	@show energy_sog, abs_error, relative_error

	# abs_error = abs(energy_ewald - energy_sog)
	# relative_error = abs(energy_ewald - energy_sog) / abs(energy_ewald)

	# @show energy_sog, abs_error, relative_error

	df = DataFrame(n_atoms = n_atoms, E_exact = energy_ewald, E_fssog = energy_sog, abs_error = abs_error, relative_error = relative_error)
	CSV.write("data/Acc_T6_cube_total_6.csv", df, append = true)
end
