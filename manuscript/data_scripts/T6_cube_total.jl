using ExTinyMD, FastSpecSoG, EwaldSummations, Plots, Random
using CSV, DataFrames, LaTeXStrings
Random.seed!(1234)

n_atoms0 = 8
L0 = 20.0
N_real0 = (32, 32, 32)

w = (12, 12, 12)
β = 5.0 .* w
extra_pad_ratio = 2
cheb_order = 10
preset = 2
M_mid = 4
N_grid = (8, 8, 6)
Q = 16
r_c = 9.99

α0 = 0.61

for ratio in [2^i for i in 0:0.5:2]

    n_atoms = Int(ceil(n_atoms0 * ratio^3))
    if isodd(n_atoms)
        n_atoms += 1
    end
    L = L0 * ratio

    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    @info "Ewald2D, n_atoms = $n_atoms"
    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 6.0, α0 / ratio, (L, L, L), ϵ = 1.0)
    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    @time energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)


    N_real = Int.(ceil.(ratio .* N_real0))

    @info "FSSoG, n_atoms = $n_atoms"
    fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Q, 0.5, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Q, 64; preset = preset, ϵ = 1.0)
    fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
    @time energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

    abs_error = abs(energy_ewald - energy_sog)
    relative_error = abs(energy_ewald - energy_sog) / abs(energy_ewald)

    @show n_atoms, energy_ewald, energy_sog, abs_error, relative_error

    df = DataFrame(n_atoms = n_atoms, E_exact = energy_ewald, E_fssog = energy_sog, abs_error = abs_error, relative_error = relative_error)
    CSV.write("data/Acc_T6_cube_total.csv", df, append=true)
end