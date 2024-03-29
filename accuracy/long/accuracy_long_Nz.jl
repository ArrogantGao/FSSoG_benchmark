using FastSpecSoG, CSV, DataFrames
using Plots
using Random
Random.seed!(1234)

n_atoms = 100

L = (100.0, 100.0, 100.0)

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

# using the FFCT

uspara = USeriesPara(2)

errors = []

M_mid_array = [1:10...]
N_z_array = [3:3:30...]


data_file = "data/accuracy_long_Nz.csv"

for i in 1:length(M_mid_array)
    M_mid = M_mid_array[i]
    E_direct_k = long_energy_us_k(qs, poses, 30, L, uspara, M_mid + 1, length(uspara.sw))

    for N_z in N_z_array
        N_grid = (32, 32, N_z)    

        k_x, k_y, r_z, H_r, H_c, phase_x, phase_y,  sort_z, z = long_paras_gen(L, N_grid, n_atoms)

        E_FFCT_loop_k = energy_long_loop_k(qs, poses, L, M_mid, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, uspara)

        abs_error = abs(E_FFCT_loop_k - E_direct_k)
        relative_error = abs(E_FFCT_loop_k - E_direct_k) / abs(E_direct_k)

        @show N_grid, M_mid, abs_error, relative_error

        df = DataFrame(N_z = N_z, M_mid = M_mid, E_exact = E_direct_k, E_FFCT = E_FFCT_loop_k, abs_error = abs_error, relative_error = relative_error)

        CSV.write(data_file, df, append=true)
    end
end

df0 = CSV.read("data/accuracy_long_Nz.csv", DataFrame)
M_mid_array = unique(df0.M_mid)
plot(xlabel = "N_z", ylabel = "log10(Relative error)", legend = :topright, dpi = 500)
for i in 1:2:length(M_mid_array)
    M_mid = M_mid_array[i]
    df = filter(row -> row.M_mid == M_mid, df0)
    plot!(df.N_z, log10.(df.relative_error), label = "M_mid = $M_mid", marker = :circle)
end
savefig("figs/accuracy_long_Nz.png")