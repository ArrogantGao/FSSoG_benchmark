using CSV, DataFrames, CairoMakie, LaTeXStrings

df_xy = CSV.read("data/Acc_T4_cube_Nxy.csv", DataFrame)
df_z = CSV.read("data/Acc_T4_cube_Rz.csv", DataFrame)

M_mid = unique(df_xy.M_mid)
N_xy = unique(df_xy.Nxy)
R_z = unique(df_z.Rz)

marker = [:circle, :diamond, :star5, :utriangle, :hexagon, :xcross]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axl = Axis(ga[1, 1], xlabel = L"N_{xy}", ylabel = L"\mathcal{E}_r", yscale = log10)
axr = Axis(gb[1, 1], xlabel = L"R_z", yscale = log10)

xlims!(axl, (0, 17))
xlims!(axr, (0, 65))
ylims!(axl, (1e-16, 1))
ylims!(axr, (1e-16, 1))

for i in 1:4
    M = M_mid[i]
    mask = df_xy.M_mid .== M
    scatter!(axl, N_xy, df_xy.error_rel[mask], markersize = 10, label = L"M_{\text{mid}} = %$M", marker = marker[i])
end

for i in 1:4
    M = M_mid[i]
    mask = df_z.M_mid .== M
    scatter!(axr, R_z, df_z.error_rel[mask], markersize = 10, label = L"M_{\text{mid}} = %$M", marker = marker[i])
end

axislegend(axl, position = :rt)

text!(axl, (3, 1e-13), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (13, 1e-13), text = "(b)", fontsize = 30, align = (:right, :baseline),)

save("figs/long_Nxyz.pdf", f)
save("figs/long_Nxyz.png", f)
f