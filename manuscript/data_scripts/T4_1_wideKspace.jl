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
uspara_array = [USeriesPara(i...) for i in 1:6]
uspara = uspara_array[6]
M_mid_array = [1:2:11...]
M_total = length(uspara.sw)
E_kspace=[]
# E_direct_total = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_total) + long_energy_us_0(qs, poses, L, uspara, 1, M_total)
E_direct_total = 293.04310539443907

E_direct_ks_l = zeros(Float64, M_total)
accuracy = 1e-16
Threads.@threads for l in 1:M_total
    s, w = uspara.sw[l]
    km = sqrt(-4 * log(accuracy) / s^2)
    cutoff = ceil(Int, km * maximum(L) / 2π) + 1
    E_direct_ks_l[l] = FastSpecSoG.long_energy_sw_k(qs, poses, cutoff, L, s, w)
    @show l, E_direct_ks_l[l]
end

E_direct_ks = [sum(E_direct_ks_l[i + 1:end]) for i in M_mid_array]
@show E_direct_ks

E_sums = zeros(Float64, (M_total, length(K_max_array)))

for j in 1:length(K_max_array)
    Threads.@threads for i in 1:M_total
        kk = K_max_array[j]
        s, w = uspara.sw[i]
        E_sums[i, j] = FastSpecSoG.long_energy_sw_k(qs, poses, kk, L, s, w)
        @show i, j, E_sums[i, j]
    end
end

@show E_sums

E_sum_k = [[sum(E_sums[i+1:end, j]) for i in M_mid_array] for j in 1:length(K_max_array)]
@show E_sum_k

for j in 1:length(M_mid_array)
    M_mid = M_mid_array[j]
    for i in 1:length(K_max_array)
        kk = K_max_array[i]

        error = E_sum_k[i][j] - E_direct_ks[j]
        error_rel = error/abs(E_direct_total)

        @show kk, M_mid, error, error_rel
        df = DataFrame(k_max = kk, energy_exact = E_direct_ks[j], abs_error = error, relative_error = error_rel, M_mid = M_mid)
        CSV.write("data/Acc_T4_1_wideKspace.csv", df, append = true)
    end
end

