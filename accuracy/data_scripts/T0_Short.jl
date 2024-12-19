using FastSpecSoG, BenchmarkTools, ExTinyMD
using Plots, CSV, DataFrames
using Random
Random.seed!(123)

@show Threads.nthreads(), nprocs()

U(q1, q2, r_sq, s, w) = w * exp(- r_sq / s^2) * q1 * q2
function U_total(qs, neighbors, s, w)
    U_total = 0.0
    for (i, j, r_sq) in neighbors
        U_total += U(qs[i], qs[j], r_sq, s, w)
    end
    for q in qs
        U_total += w * q^2 / 2
    end
    return U_total
end

CSV.write("data/Acc_T0_short.csv", DataFrame(preset = Int[], M = Int[], total = Float64[], near_sw = Float64[], error_rel = Float64[]))

n_atoms = 1000
L = (20.0, 20.0, 20.0)

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
r_c = 9.9999999

neighbors = []
boundary = ExTinyMD.Q2dBoundary(L...)
for i in 1:n_atoms - 1
    for j in i + 1:n_atoms
        if i != j
            pos_1, pos_2, r_sq = position_check3D(poses[i], poses[j], boundary, r_c)
            if r_sq != 0.0
                push!(neighbors, (i, j, r_sq))
            end
        end
    end
end

M_max = 250

accuracy = 1e-16

# energy_totals = []

# for preset in 1:6
#     uspara0 = USeriesPara(preset)
#     M_t = length(uspara0.sw)
#     energy_total = long_energy_us_k(qs, poses, 1e-16, L, uspara0, 1, M_t) + long_energy_us_0(qs, poses, L, uspara0, 1, M_t)
#     push!(energy_totals, energy_total)
#     @show preset, energy_total
# end
# @show energy_totals

energy_totals = [96.63438269437341, 124.4380320000621, 164.66320580651364, 183.8203574798731, 229.872932462265, 293.043105394439]

total_contribution = zeros(6, M_max)
near_sw = zeros(6, M_max)
for preset in 1:6
    b, σ, ω, M0 = FastSpecSoG.preset_parameters[preset]
    uspara = USeriesPara(b, σ, ω, M_max)
    Threads.@threads for i in 1:8:M_max
        s, w = uspara.sw[i]
        km = sqrt(-4 * log(accuracy) / s^2)
        cutoff = ceil(Int, km * maximum(L) / 2π) + 1
        total_contribution[preset, i] = long_energy_sw_k(qs, poses, cutoff, L, s, w) + long_energy_sw_0(qs, poses, L, s, w)
        near_sw[preset, i] = U_total(qs, neighbors, s, w)
        error_rel = abs(total_contribution[preset, i] - near_sw[preset, i]) / abs(energy_totals[preset])
        @show preset, i, total_contribution[preset, i], near_sw[preset, i], error_rel
    end
end

for preset in 1:6
    for i in 1:8:M_max
        df = DataFrame(preset = preset, M = i, total = total_contribution[preset, i], near_sw = near_sw[preset, i], error_rel = abs(total_contribution[preset, i] - near_sw[preset, i]) / abs(energy_totals[preset]))
        CSV.write("data/Acc_T0_short.csv", df, append = true)
    end
end