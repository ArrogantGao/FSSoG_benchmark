using FastSpecSoG, CSV, DataFrames, EwaldSummations, ExTinyMD, ForwardDiff, JLD2
using Base.Threads
using Random
Random.seed!(1234)

function USeries_energy(q_1::T, q_2::T, x_1::Point{3, T}, x_2::Point{3, T}, uspara::USeriesPara{T}) where{T}
    Δr = sqrt(dist2(x_1, x_2))
    if Δr ≥ 10.0
        return q_1 * q_2 * (1/Δr - U_series(Δr, uspara)) / 4π
    else
        return zero(T)
    end
end

@load "reference/cube/n_100_energy.jld2" n_atoms L atoms info energy_ewald

poses = [info.particle_info[i].position for i in 1:n_atoms]
qs = [atoms[i].charge for i in 1:n_atoms]

for preset in 1:6
    for (iM, M) in enumerate([1:15:300...])
        b, σ, ω, M0 = FastSpecSoG.preset_parameters[preset]
        uspara = USeriesPara(b, σ, ω, M)
        energy_total = 0.0
        for i in 1:n_atoms
            ei = Atomic{Float64}(0.0)
            x_1 = poses[i]
            @threads for j in 1:n_atoms
                t = 0.0
                for mx in -40:40, my in -40:40
                    x_2 = poses[j] + Point(20.0 * mx, 20.0 * my, 0.0)
                    t += USeries_energy(qs[i], qs[j], x_1, x_2, uspara)
                end
                atomic_add!(ei, t)
            end
            energy_total += ei[]
        end
        e = abs(energy_total) / abs(energy_ewald)
        @show preset, M, e
        df = DataFrame(preset = preset, M = M, error = e)
        CSV.write("data/Acc_T0_energy.csv", df, append = true)
    end
end