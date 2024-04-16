using CSV, DataFrames, CairoMakie, LaTeXStrings, LsqFit, FastSpecSoG

df_energy = CSV.read("data/Acc_T0_energy.csv", DataFrame)
df_force = CSV.read("data/Acc_T0_force.csv", DataFrame)

Ms = unique(df_force.M)
Mse = unique(df_energy.M)
presets = unique(df_force.preset)

marker = [:circle, :diamond, :star5, :utriangle, :hexagon, :xcross]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axl = Axis(ga[1, 1], xlabel = L"M", ylabel = L"\mathcal{E}_r", yscale = log10)
axr = Axis(gb[1, 1], xlabel = L"M", yscale = log10)

xlims!(axl, (0, 300))
xlims!(axr, (0, 120))

ylims!(axl, (1e-16, 1))
ylims!(axr, (1e-16, 1))

for preset in presets
    mask_energy = df_energy.preset .== preset
    mask_force = df_force.preset .== preset
    b = round(FastSpecSoG.preset_parameters[preset][1], digits = 2)
    scatter!(axl, Mse, df_energy.error[mask_energy], markersize = 10, label = L"b = %$b", marker = marker[preset])
    scatter!(axr, Ms, df_force.error[mask_force], markersize = 10, label = L"b = %$b", marker = marker[preset])
end

axislegend(axl, position = :lb, nbanks = 2)

text!(axl, (50, 1e-10), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (25, 1e-10), text = "(b)", fontsize = 30, align = (:right, :baseline),)

f

save("figs/cutoff.pdf", f)
save("figs/cutoff.png", f)