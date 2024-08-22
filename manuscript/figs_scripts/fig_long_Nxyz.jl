using CSV, DataFrames, CairoMakie, LaTeXStrings, LsqFit, FastSpecSoG
include("setting.jl")

df_xy = CSV.read("data/Acc_T4_cube_Nxy.csv", DataFrame)
df_z = CSV.read("data/Acc_T4_cube_Rz.csv", DataFrame)

M_mid = unique(df_xy.M_mid)
N_xy = unique(df_xy.Nxy)
R_z = unique(df_z.Rz)

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axl = Axis(ga[1, 1], xlabel = L"I^{\mathcal{M}_\text{l}}", ylabel = L"\mathcal{E}_r", yscale = log10)
axr = Axis(gb[1, 1], xlabel = L"P", yscale = log10)

xlims!(axl, (0, 35))
xlims!(axr, (0, 65))
ylims!(axl, (1e-16, 1))
ylims!(axr, (1e-16, 1))

eta_array=[0.11,0.15,0.19,0.26]

for i in 1:4
    M = M_mid[i]
    mask = df_xy.M_mid .== M
    eta = eta_array[i]
    scatter!(axl, N_xy, df_xy.error_rel[mask], markersize = ms, label = L"\eta \approx %$eta", marker = markers[i], color = colors[i])

    s, w = USeriesPara(6).sw[M + 1]
    @. model(x, p) = log.(abs(p[1] * w * exp(- s^2 * p[2]* x^2)))

    ks = range(0, 32, length = 1000)

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
    lines!(axl, ks, exp.(g.(ks)), linestyle = :dash, linewidth = lw, color = colors[i])
end

for i in 1:4
    M = M_mid[i]
    mask = df_z.M_mid .== M
    eta = eta_array[i]
    scatter!(axr, R_z, df_z.error_rel[mask], markersize = ms, label = L"\eta \approx %$eta", marker = markers[i], color = colors[i])

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
    lines!(axr, xs, exp.(g.(xs)), linestyle = :dash, linewidth = lw, color = colors[i])
end

axislegend(axl, position = :rt)

text!(axl, (7, 1e-13), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (13, 1e-13), text = "(b)", fontsize = 30, align = (:right, :baseline),)

text!(axl, (33, 10^(-7.8)), text = L"O(e^{- C_1 (I^{\mathcal{M}_\text{l}})^2})", fontsize = 20, align = (:right, :baseline),)
text!(axr, (63, 10^(-7.8)), text = L"O(C_2^{-P} / \sqrt{P !})", fontsize = 20, align = (:right, :baseline),)

save("figs/long_Nxyz.pdf", f)
save("figs/long_Nxyz.png", f)
f