using FastSpecSoG, CSV, DataFrames
using CairoMakie, LaTeXStrings

df_cube = CSV.read("data/Acc_T4_cube_zeroth.csv", DataFrame)
df_thin = CSV.read("data/Acc_T5_thin_zeroth.csv", DataFrame)

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axl = Axis(ga[1, 1], xlabel = L"Q_0", ylabel = L"\mathcal{E}_r", yscale = log10)
axr = Axis(gb[1, 1], xlabel = L"Q_0", yscale = log10)

Q_l = unique(df_cube.Q_0)
Q_r = unique(df_thin.Q_0)


M_mids = unique(df_cube.M_mid)
gammas = unique(df_thin.gamma)

marker = [:circle, :diamond, :star5, :utriangle, :hexagon, :xcross]

Is = [1, 3, 5]

for i0 in 1:3
    i = Is[i0]
    M_mid = M_mids[i]
    mask = df_cube.M_mid .== M_mid
    scatter!(axl, Q_l, abs.(df_cube.error_rel[mask] .+ 1e-16), markersize = 10, label = L"M_{\text{mid}} = %$M_mid", marker = marker[i0])
end

for i in 1:length(gammas)
    gamma = gammas[i]
    mask = df_thin.gamma .== gamma
    scatter!(axr, Q_r, abs.(df_thin.error_rel[mask] .+ 1e-16), markersize = 10, label = L"\gamma = %$gamma", marker = marker[i])
end

axislegend(axl, position = :rt)
axislegend(axr, position = :rt)

text!(axl, (12, 1e-14), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (4, 1e-14), text = "(b)", fontsize = 30, align = (:right, :baseline),)

save("figs/zeroth.pdf", f)
save("figs/zeroth.png", f)

f