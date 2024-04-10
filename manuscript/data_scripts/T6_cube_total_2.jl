using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

L0 = 20.0

N_real0 = (32, 32, 32)
w = (16, 16, 16)
β = 5.0 .* w
extra_pad_ratio = 2
cheb_order = 10
preset = 2
M_mid = 4
N_grid = (8, 8, 6)
Q = 16
r_c = 9.99

for data in ["n_1000.jld2", "n_3164.jld2", "n_10000.jld2", "n_31624.jld2", "n_100000.jld2"]

    path = "./reference/cube/$data" 
	@load path n_atoms L atoms info energy_ewald

	@show n_atoms, L, energy_ewald

	boundary = ExTinyMD.Q2dBoundary(L, L, L)

    ratio = (n_atoms / 1000)^(1/3)
	N_real = Int.(ceil.(ratio .* N_real0))

	@info "FSSoG, n_atoms = $n_atoms"
	fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Q, 0.5, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Q, 64; preset = preset, ϵ = 1.0)
	fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
	@time energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

	abs_error = abs(energy_ewald - energy_sog)
	relative_error = abs(energy_ewald - energy_sog) / abs(energy_ewald)

	@show energy_sog, abs_error, relative_error

	df = DataFrame(n_atoms = n_atoms, E_exact = energy_ewald, E_fssog = energy_sog, abs_error = abs_error, relative_error = relative_error)
	CSV.write("data/Acc_T6_cube_total_2.csv", df, append = true)
end