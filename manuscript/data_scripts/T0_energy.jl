using FastSpecSoG, CSV, DataFrames, EwaldSummations, ExTinyMD, ForwardDiff, JLD2
using OhMyThreads: tmapreduce
using Base.Threads: nthreads
using Random
Random.seed!(1234)

@inline function USeries_energy(q_1::T, q_2::T, x_1::Point{3, T}, x_2::Point{3, T}, uspara::USeriesPara{T}, r_c::T) where{T}
    Δr = sqrt(dist2(x_1, x_2))
    return Δr ≥ r_c ? q_1 * q_2 * (1/Δr - U_series(Δr, uspara)) / 4π : zero(T)
end

@load "reference/cube_L1/n_100_energy.jld2" n_atoms L atoms info energy_ewald

# CSV.write("data/Acc_T0_energy.csv", DataFrame(preset = Int[], M = Int[], error = Float64[]))

poses = [info.particle_info[i].position for i in 1:n_atoms]
qs = [atoms[i].charge for i in 1:n_atoms]

Lx = 1.0
Ly = 1.0
Lz = 1.0
rc = 0.1

for preset in 1:6
    for (iM, M) in enumerate([286:15:350...])
        b, σ0, ω, M0 = FastSpecSoG.preset_parameters[preset]
        σ = σ0 * rc
        uspara = USeriesPara(b, σ, ω, M)

        energy_total = 0.0

        for i in 1:n_atoms
            ei = tmapreduce(+, 1:n_atoms; ntasks = nthreads()) do j
                t = 0.0
                for mx in -40:40, my in -40:40
                    x_2 = poses[j] + Point(Lx * mx, Ly * my, 0.0)
                    t += USeries_energy(qs[i], qs[j], poses[i], x_2, uspara, 0.1)
                end
                t
            end

            energy_total += ei
        end
        e = abs(energy_total) / abs(energy_ewald)
        @show preset, M, e
        df = DataFrame(preset = preset, M = M, error = e)
        CSV.write("data/Acc_T0_energy.csv", df, append = true)
    end
end