using FastSpecSoG, CSV, DataFrames, EwaldSummations, ExTinyMD, ForwardDiff, JLD2
using Base.Threads
using Random
Random.seed!(1234)

function dU_series(x::T, uspara::USeriesPara{T}) where{T}
    dU_series_value = zero(T)

    for i in eachindex(uspara.sw)
        (s, w) = uspara.sw[i]
        dU_series_value += - 2 * x * w / s^2 * exp(- x^2 / s^2)
    end

    return dU_series_value
end

function USeries_force(q_1::T, q_2::T, x_1::Point{3, T}, x_2::Point{3, T}, uspara::USeriesPara{T}) where{T}
    f = r -> q_1 * q_2 * (-1/r^2 - dU_series(r, uspara)) / 4π
    Δr = sqrt(dist2(x_1, x_2))
    if Δr ≥ 10.0
        force = - f(Δr)
        return force * (x_1 - x_2) / Δr
    else
        return Point(zero(T), zero(T), zero(T))
    end
end

function error(f, fe)
    error = 0.0
    for i in eachindex(f)
        error += dist2(f[i], Point(0.0, 0.0, 0.0)) / dist2(fe[i], Point(0.0, 0.0, 0.0))
    end
    return sqrt(error / length(f))
end

@load "reference/cube/n_100_force.jld2" n_atoms L atoms info force_ewald

force_exact = [Point(f) for f in force_ewald]

poses = [info.particle_info[i].position for i in 1:n_atoms]
qs = [atoms[i].charge for i in 1:n_atoms]

for preset in 1:6
    for (iM, M) in enumerate([1:5:250...])
        b, σ, ω, M0 = FastSpecSoG.preset_parameters[preset]
        uspara = USeriesPara(b, σ, ω, M)
        f = [Point(0.0, 0.0, 0.0) for i in 1:n_atoms]
        for i in 1:n_atoms
            x_1 = poses[i]
            fijx = Atomic{Float64}(0.0)
            fijy = Atomic{Float64}(0.0)
            fijz = Atomic{Float64}(0.0)
            @threads for j in 1:n_atoms
                fj = Point(0.0, 0.0, 0.0)
                for mx in -20:20, my in -20:20
                    x_2 = poses[j] + Point(20.0 * mx, 20.0 * my, 0.0)
                    fij = USeries_force(qs[i], qs[j], x_1, x_2, uspara)
                    fj += fij
                end
                atomic_add!(fijx, fj.coo[1])
                atomic_add!(fijy, fj.coo[2])
                atomic_add!(fijz, fj.coo[3])
            end
            f[i] += Point(fijx[], fijy[], fijz[])
        end
        e = error(f, force_exact)
        df = DataFrame(preset = preset, M = M, error = e)
        CSV.write("data/Acc_T0_force.csv", df, append = true)
    end
end