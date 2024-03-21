using FastSpecSoG, CSV, DataFrames, LaTeXStrings
using Plots
using Random
Random.seed!(1234)

n_atoms = 100

L = (100.0, 100.0, 100.0)

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]
data_file = "data/accuracy_long_Mmid.csv"

for preset in 1:5

    uspara = USeriesPara(preset)

    M_mid_array = [1:length(uspara.sw) - 1...]

    for M_mid in M_mid_array
        E_direct_k = long_energy_us_k(qs, poses, 30, L, uspara, M_mid, M_mid + 1) 
        E_direct_0 = long_energy_us_0(qs, poses, L, uspara, M_mid, M_mid + 1)
        df = DataFrame(preset = preset, M_mid = M_mid, E_exact = E_direct_k, E_0 = E_direct_0)
        CSV.write(data_file, df, append=true)
        @show preset, M_mid, E_direct_k, E_direct_0
    end
end

df = CSV.read("data/accuracy_long_Mmid.csv", DataFrame)
plot(xlabel = L"M_{mid}", ylabel = "Error", title = L"Contribution of long-range part to $M_{mid}$", dpi = 500, xlim = [0, 20.5], ylim = [-15, 1.0])
for i in 1:5
    dfi = filter(row -> row.preset == i, df)
    plot!(dfi.M_mid, log10.(abs.(dfi.E_exact)), label = "preset = $(i)", marker = :circle)
end
plot!(legend = :bottomright)
savefig("figs/accuracy_long_Mmid.png")