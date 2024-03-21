using FastSpecSoG, Plots, CSV, DataFrames, Random, ExTinyMD, LaTeXStrings
Random.seed!(1234)

N_grids = [(16, 16, 16), (32, 32, 32), (64, 64, 64), (128, 128, 128), (256, 256, 256)]

n_atoms = 100
L = (100.0, 100.0, 100.0)
w = (16, 16, 16)
β = 5.0 .* w
extra_pad_ratios = [2, 3, 4, 8]
cheb_order = 10
uspara = USeriesPara(2)
M_mid = 6

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

E_direct_ks = long_energy_us_k(qs, poses, 50, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid)

errors = [[] for j in 1:length(extra_pad_ratios)]
E_mid = [[] for j in 1:length(extra_pad_ratios)]

for i in 1:length(N_grids)
    for j in 1:length(extra_pad_ratios)
        N_real = N_grids[i]
        extra_pad_ratio = extra_pad_ratios[j]
        
        gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
        E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

        error = abs(E_3DFFT - E_direct_ks)
        push!(errors[j], error)
        push!(E_mid[j], E_3DFFT)

        @show N_real, extra_pad_ratio, error
    end
end

Ns = [N_grid[1] for N_grid in N_grids]

df = DataFrame(extra_pad_ratio = repeat(extra_pad_ratios, inner = length(Ns)), N_real = repeat(Ns, outer = length(extra_pad_ratios)), error = vcat(errors...))
CSV.write("data/accuracy_mid_Ngrid.csv", df)

plot(xlabel = L"log2(N_{real})", ylabel = L"Error", title = L"$N_{real}$ $with$ $w = (16, 16, 16)$", dpi = 500)
for i in 1:length(extra_pad_ratios)
    plot!(log2.([N_grid[1] for N_grid in N_grids]), log10.(abs.(errors[i])), label = "λ = $(extra_pad_ratios[i])", marker = :circle, xticks = [4, 5, 6, 7, 8])
end

savefig("figs/accuracy_mid_Ngrid.png")