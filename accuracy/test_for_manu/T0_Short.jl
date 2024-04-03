using FastSpecSoG, BenchmarkTools, ExTinyMD, Plots, LaTeXStrings
using Plots, CSV, DataFrames
using Random
Random.seed!(123)

n_atoms = 1000
L = 40.0
boundary = ExTinyMD.Q2dBoundary(L, L, L)

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms÷2
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
end

for i in n_atoms÷2 + 1 : n_atoms
    push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)
position = [info.particle_info[i].position.coo for i in 1:n_atoms]
charge = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]

# charge = rand(n_atoms)
# charge .-= sum(charge) ./ n_atoms

energy_short_M = [[] for i in 1:6]
cheb_orders = [1:3:40...]

# fix r_c = 10.0
r_c = 10.0

neighbor = CellList3D(info, r_c, boundary, 1)

for preset in 1:6
    b, σ, ω, M0 = FastSpecSoG.preset_parameters[preset]
    for M in 1:250
        interaction_short = FSSoG_naive((L, L, L), n_atoms, r_c, 3.0, b, σ, ω, M)
        Es = short_energy_naive(interaction_short, neighbor, position, charge)
        push!(energy_short_M[preset], Es)
    end
end

contribution = [[] for i in 1:6]
for i in 1:6
    for j in 1:length(energy_short_M[i]) - 1
        push!(contribution[i], abs(energy_short_M[i][end] - energy_short_M[i][j]))
        df = DataFrame(uspara = i, error = abs(energy_short_M[i][end] - energy_short_M[i][j]) /abs(energy_short_M[i][end]))
        CSV.write("data_for_manu/Acc_T0_Short.csv", df, append = true)
    end
end

