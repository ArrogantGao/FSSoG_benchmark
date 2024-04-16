using CSV, DataFrames, CairoMakie, LaTeXStrings, LsqFit, FastSpecSoG

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

    s, w = USeriesPara(6).sw[M + 1]
    @. model(x, p) = log.(abs(p[1] * w * exp(- s^2 * p[2]* x^2)))

    ks = range(0, 16, length = 1000)

    raw_y_data = df_xy.error_rel[mask]
    raw_x_data = N_xy
    n0 = 4
    n = findfirst(x -> x < 1e-14, raw_y_data)
    x_data = raw_x_data[n0:n]
    y_data = raw_y_data[n0:n]
    p0 = [1.0, 0.01]
    fit = curve_fit(model, x_data, log.(y_data), p0)

    @info "(a) $i $(fit.param[1] * w) exp(-x^2 $(s^2 * fit.param[2]))"

    g = x -> model(x, fit.param)
    lines!(axl, ks, exp.(g.(ks)), linestyle = :dash, linewidth = 0.7)
end

for i in 1:4
    M = M_mid[i]
    mask = df_z.M_mid .== M
    scatter!(axr, R_z, df_z.error_rel[mask], markersize = 10, label = L"M_{\text{mid}} = %$M", marker = marker[i])

    @. model(x, p) = log.(abs(p[1] / (p[2]^x * sqrt(factorial(big(x)))) ))

    xs = [1:65...]

    raw_y_data = df_z.error_rel[mask]
    raw_x_data = R_z
    n0 = 1
    n = findfirst(x -> x < 1e-14, raw_y_data)
    x_data = raw_x_data[n0:n]
    y_data = raw_y_data[n0:n]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, x_data, log.(y_data), p0)

    @info "(b) $i $(fit.param[1]) / ( $(fit.param[2])^x * sqrt(x!))"

    g = x -> model(x, fit.param)
    lines!(axr, xs, exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

axislegend(axl, position = :rt)

text!(axl, (3, 1e-13), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (13, 1e-13), text = "(b)", fontsize = 30, align = (:right, :baseline),)

text!(axl, (15, 10^(-7.8)), text = L"O(e^{- C_1 N_{xy}^2})", fontsize = 20, align = (:right, :baseline),)
text!(axr, (63, 10^(-7.8)), text = L"O(C_2^{-R_z} / \sqrt{R_z !})", fontsize = 20, align = (:right, :baseline),)

save("figs/long_Nxyz.pdf", f)
save("figs/long_Nxyz.png", f)
f