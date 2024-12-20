using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

L0 = 20.0
N_real0 = (48, 48, 48)

w = (12, 12, 12)
β = 5.0 .* w
extra_pad_ratio_intial = 10
cheb_order = 16
preset = 6
uspara = USeriesPara(preset)
M_mid_initial = 15
eta = uspara.sw[M_mid_initial][1] / L0 + 0.0001
N_grid = (8, 8, 20)
Qs = 32
Ql = 32
r_c = 9.99
Q_0 = 32
Rz_0 = 32

for data in ["n_1000.jld2", "n_3164.jld2", "n_10000.jld2", "n_31624.jld2", "n_100000.jld2"]

    path = "./reference/cube/$data" 
	@load path n_atoms L atoms info energy_ewald

	@show n_atoms, L, energy_ewald

	boundary = ExTinyMD.Q2dBoundary(L, L, L)

    ratio = (n_atoms / 1000)^(1/3)
	N_real = Int.(ceil.(ratio .* N_real0))
	extra_pad_ratio = Int(ceil(extra_pad_ratio_intial * ratio))

	M_mid = proper_M(eta, L, uspara)
	@show M_mid

	@info "FSSoG, n_atoms = $n_atoms"
	fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Qs, 0.5, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Ql, Rz_0, Q_0; preset = preset, ϵ = 1.0)
	fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
	@time energy_fssog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

	abs_error = abs(energy_ewald - energy_fssog)
	relative_error = abs(energy_ewald - energy_fssog) / abs(energy_ewald)

	@show energy_fssog, abs_error, relative_error

    @save path n_atoms L atoms info energy_ewald energy_fssog

	# df = DataFrame(n_atoms = n_atoms, E_exact = energy_ewald, E_fssog = energy_sog, abs_error = abs_error, relative_error = relative_error)
	# CSV.write("data/Acc_T6_cube_total_6.csv", df, append = true)
end
