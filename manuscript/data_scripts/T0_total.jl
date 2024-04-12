using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

L0 = 20.0
N_real0 = (36, 36, 36)

w = (8, 8, 8)
β = 5.0 .* w
extra_pad_ratio_intial = 8
cheb_order = 10
M_mid_initial = 15
uspara = USeriesPara(6)
eta = uspara.sw[M_mid_initial][1] / L0 + 0.0001
N_grid = (2, 2, 15)
Qs = 28
Ql = 12
r_c = 9.99
Q_0 = 32
Rz_0 = 16

for preset in 1:6
    uspara = USeriesPara(preset)
    path = "./reference/cube/n_3164.jld2" 
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
	@time energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

	abs_error = abs(energy_ewald - energy_sog)
	relative_error = abs(energy_ewald - energy_sog) / abs(energy_ewald)

	@show preset, energy_sog, abs_error, relative_error

	df = DataFrame(preset = preset, E_exact = energy_ewald, E_fssog = energy_sog, abs_error = abs_error, relative_error = relative_error)
	CSV.write("data/Acc_T0_total.csv", df, append = true)
end
