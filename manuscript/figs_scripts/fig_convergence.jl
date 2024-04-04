using CSV, DataFrames, CairoMakie, LaTeXStrings, FastSpecSoG

df_k = CSV.read("data/Acc_T1_Kspace.csv", DataFrame)
df_r = CSV.read("data/Acc_T0_short.csv", DataFrame)

L = 20.0
km = unique(df_k.k_max) .* π ./ L
Ms = unique(df_r.M_mid)

marker = [:circle, :diamond, :star5, :utriangle, :hexagon, :xcross]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axr = Axis(gb[1, 1], xlabel = L"k_{\text{max}}", yscale = log10)
axl = Axis(ga[1, 1], xlabel = L"M", ylabel = L"\mathcal{E}_r", yscale = log10)

xlims!(axr, (0, 1.1 * π))
xlims!(axl, (5, 250))
ylims!(axr, (1e-16, 1))
ylims!(axl, (1e-16, 1))

for i in 1:6
    mask_1 = df_k.preset .== i
    mask_2 = df_r.preset .== i
    b = round(FastSpecSoG.preset_parameters[i][1], digits = 2)
    scatter!(axr, km, df_k.error_rel[mask_1] .+ 3e-16, markersize = 10, label = L"b = %$b", marker = marker[i])
    scatter!(axl, Ms, df_r.error_rel[mask_2], markersize = 10, label = L"b = %$b", marker = marker[i])
end

axislegend(axr, position = :rt)

text!(axl, (35, 1e-14), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (0.5, 1e-14), text = "(b)", fontsize = 30, align = (:right, :baseline),)

save("figs/convergence.pdf", f)



f