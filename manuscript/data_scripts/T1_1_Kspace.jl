using FastSpecSoG, CSV, DataFrames
using Random
Random.seed!(123)

@show Threads.nthreads(), nprocs()

n_atoms = 1000
L = (20.0, 20.0, 20.0)

#qs = rand(n_atoms)
#qs .-= sum(qs) ./ n_atoms
qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
K_max_array = [1:20...]
M_mid = 0

for preset in 1:6
    uspara = USeriesPara(preset)
    M_total = length(uspara.sw)
    E_kspace=[]
    E_direct_total_k = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_total)
    E_direct_total_0 = long_energy_us_0(qs, poses, L, uspara, 1, M_total)
    E_direct_total = E_direct_total_k + E_direct_total_0
    @show E_direct_total, E_direct_total_k, E_direct_total_0

    E_sums = zeros(Float64, (M_total, length(K_max_array)))

    Threads.@threads for (i, j) in [(i, j) for i in 1:M_total, j in 1:length(K_max_array)]
        kk = K_max_array[j]
        s, w = uspara.sw[i]
        E_sums[i, j] = FastSpecSoG.long_energy_sw_k(qs, poses, kk, L, s, w)
        @show i, j, E_sums[i, j]
    end

    E_sum_k = [sum(E_sums[1:end, j]) for j in 1:length(K_max_array)]
    @show E_sum_k

    for i in 1:length(K_max_array)
        kk = K_max_array[i]

        error = abs(E_sum_k[i] - E_direct_total_k)
        error_rel = error/abs(E_direct_total)

        @show kk, preset, error, error_rel
        df = DataFrame(preset = preset, k_max = kk, energy_exact = E_direct_total, energy_k = E_sum_k[i] + E_direct_total_0, abs_error = error, relative_error = error_rel)
        CSV.write("data/Acc_T1_Kspace.csv", df, append = true)
    end
end

