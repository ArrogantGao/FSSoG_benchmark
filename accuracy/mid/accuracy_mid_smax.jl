using FastSpecSoG, Plots, CSV, DataFrames, LaTeXStrings

using Random
Random.seed!(1234)

n_atoms = 100
L = (100.0, 100.0, 100.0)
N_real = (128, 128, 128)
w = (16, 16, 16)
β = 5.0 .* w
extra_pad_ratios = [2, 3, 4, 8]
cheb_order = 10
uspara = USeriesPara(2)
M_mids = [1:20...]

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

errors = [[] for i in 1:length(extra_pad_ratios)]

E_direct_ks = [long_energy_us_k(qs, poses, 50, L, uspara, 1, M_mid) + long_energy_us_0(qs, poses, L, uspara, 1, M_mid) for M_mid in M_mids]

for i in 1:length(extra_pad_ratios)
    for j in 1:length(M_mids)
        extra_pad_ratio = extra_pad_ratios[i]
        M_mid = M_mids[j]

        gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
        E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

        error = abs(E_3DFFT - E_direct_ks[j])
        
        @show extra_pad_ratio, M_mid, error

        push!(errors[i], error)
    end
end

df = DataFrame(extra_pad_ratio = repeat(extra_pad_ratios, inner = length(M_mids)), M_mid = repeat(M_mids, outer = length(extra_pad_ratios)), error = vcat(errors...))
CSV.write("data/accuracy_mid_smax.csv", df)

s_max = [uspara.sw[M_mid][1] for M_mid in M_mids]
plot(xlabel = "s_max", ylabel = "Error", title = "Error vs s_max", dpi = 500)
for i in 1:length(extra_pad_ratios)
    plot!(s_max, abs.(errors[i]), xscale = :log10, yscale = :log10, label = "extra_pad_ratio = $(extra_pad_ratios[i])", marker = :circle)
end
savefig("figs/accuracy_mid_smax.png")

plot(xlabel = L"M_{mid}", ylabel = "Error", title = L"Accuracy of mid-range part to $M_{mid}$", dpi = 500, xlim = [0, 20.5])
for i in 1:length(extra_pad_ratio)
    plot!(M_mids, errors[i], label = "λ = $(extra_pad_ratio[i])", marker = :circle, yscale = :log10)
end
plot!(legend = :bottomright)
savefig("figs/accuracy_mid_Mmid.png")