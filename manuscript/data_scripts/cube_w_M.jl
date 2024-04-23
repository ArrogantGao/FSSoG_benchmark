using FastSpecSoG, CSV, DataFrames
using Random
Random.seed!(123)

@show Threads.nthreads()

n_atoms = 100
L = (10.0, 10.0, 10.0)

qs = [(-1.0)^i for i in 1:n_atoms]
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
uspara = USeriesPara(6, r_c = 1.79)
M_mid_array = [3, 9, 15, 17, 19, 21]
w_array = [(i, i, i) for i in 1:16]
N_real = (94, 94, 94)
data_file = "data/cube_w_M.csv"

for i in 1:length(M_mid_array)
    M_mid = M_mid_array[i]
    E_direct_k = long_energy_us_k(qs, poses, 1e-16, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)

    for w in w_array
        β = 5.0 .* w
        extra_pad_ratio = 8
        cheb_order = 16
        gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
        E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)
        abs_error = abs(E_3DFFT - E_direct_k)
        rel_error = abs_error / abs(E_direct_k)

        @show w[1], M_mid, abs_error, rel_error, E_3DFFT, E_direct_k

        df = DataFrame(w = w[1], M_mid = M_mid, abs_error = abs_error, rel_error = rel_error)

        CSV.write(data_file, df, append=true)
    end
end