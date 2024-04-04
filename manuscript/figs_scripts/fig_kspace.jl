using CSV, DataFrames, CairoMakie, LaTeXStrings, FastSpecSoG

df_k = CSV.read("data/Acc_T1_Kspace.csv", DataFrame)

L = 20.0
km = unique(df_k.k_max) .* π ./ L

marker = [:circle, :diamond, :star5, :utriangle, :hexagon, :xcross]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (400, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()

axl = Axis(ga[1, 1], xlabel = L"k_{\text{max}}", ylabel = L"\mathcal{E}_r", yscale = log10)
# axr = Axis(gb[1, 1], xlabel = L"R_z", yscale = log10)

xlims!(axl, (0, 1.1 * π))
# xlims!(axr, (0, 65))
ylims!(axl, (1e-16, 1))
# ylims!(axr, (1e-16, 1))

for i in 1:6
    mask = df_k.preset .== i
    b = round(FastSpecSoG.preset_parameters[i][1], digits = 2)
    scatter!(axl, km, df_k.error_rel[mask] .+ 3e-16, markersize = 10, label = L"b = %$b", marker = marker[i])
end

axislegend(axl, position = :rt)


save("figs/kspace.pdf", f)

f