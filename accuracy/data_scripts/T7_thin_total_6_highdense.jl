using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

CSV.write("data/Acc_T7_thin_total_6_highdense.csv", DataFrame(n_atoms = Int[], abs_error = Float64[], relative_error = Float64[]))

N_real0 = (130, 130)
R_z = 8
w = (9, 9)
β = 7.5 .* w
cheb_order = 8
preset = 6
Q = 32
Q_0 = 16
Rz_0 = 12
Taylor_Q = 6
r_c = 1.3

#for data in ["n_1000.jld2", "n_3164.jld2", "n_10000.jld2", "n_31624.jld2", "n_100000.jld2"]
for data in ["n_100000.jld2"]

    path = "./reference/thin/high_rho/$data" 
	@load path n_atoms Lx Lz atoms info energy_ewald 

	energy_exact = energy_ewald

	@show n_atoms, Lx, Lz, energy_ewald

	boundary = ExTinyMD.Q2dBoundary(Lx, Lx, Lz)

    ratio = 1.0
	N_real = Int.(ceil.(ratio .* N_real0))

	@info "FSSoG, n_atoms = $n_atoms"
	fssog_interaction = FSSoGThinInteraction((Lx, Lx, Lz), n_atoms, r_c, Q, 0.01, N_real, R_z, w, β, cheb_order, Taylor_Q, Rz_0, Q_0; preset = preset, ϵ = 1.0)
	fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
	@time energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

	abs_error = abs(energy_exact - energy_sog)
	relative_error = abs(energy_exact - energy_sog) / abs(energy_exact)

	@show energy_sog, abs_error, relative_error

	@time E_short = energy_short(fssog_interaction, fssog_neighbor)
	@time E_long = energy_long(fssog_interaction)

	df = DataFrame(n_atoms = n_atoms, E_exact = energy_exact, E_fssog = energy_sog, abs_error = abs_error, relative_error = relative_error)
	CSV.write("data/Acc_T7_thin_total_6_highdense.csv", df, append = true)
end
