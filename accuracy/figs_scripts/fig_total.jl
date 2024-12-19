using CSV, DataFrames, CairoMakie, LaTeXStrings, FastSpecSoG, LsqFit, Statistics

include("setting.jl")

df_cube_2 = CSV.read("data/Acc_T6_cube_total_2.csv", DataFrame)
df_cube_4 = CSV.read("data/Acc_T6_cube_total_4.csv", DataFrame)
df_cube_6 = CSV.read("data/Acc_T6_cube_total_6.csv", DataFrame)
df_thin_2 = CSV.read("data/Acc_T6_thin_total_2.csv", DataFrame)
df_thin_4 = CSV.read("data/Acc_T6_thin_total_4.csv", DataFrame)
df_thin_6 = CSV.read("data/Acc_T6_thin_total_6.csv", DataFrame)

df_t_cube = CSV.read("data/benchmarks_cube.csv", DataFrame)
df_t_thin = CSV.read("data/benchmarks_thin.csv", DataFrame)

df_cube = [df_cube_2, df_cube_4, df_cube_6]
df_thin = [df_thin_2, df_thin_4, df_thin_6]

n_atoms = [1000, 3164, 10000, 31624, 100000]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 700), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

gta = f[2, 1] = GridLayout()
gtb = f[2, 2] = GridLayout()

axr = Axis(gb[1, 1], xlabel = L"N", yscale = log10, xscale = log10, xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"]))
axl = Axis(ga[1, 1], xlabel = L"N", ylabel = L"\mathcal{E}_r", yscale = log10, xscale = log10, xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"]))

axtl = Axis(gta[1, 1], xlabel = L"N", ylabel = L"\text{Time (s)}", yscale = log10, xscale = log10, yticks = ([1e-2, 0.1, 1, 10], ["10⁻²", "10⁻¹", "10⁰", "10¹"]), xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"]))
axtr = Axis(gtb[1, 1], xlabel = L"N", yscale = log10, xscale = log10, yticks = ([1e-2, 0.1, 1, 10], ["10⁻²", "10⁻¹", "10⁰", "10¹"]), xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"]))

ylims!(axr, (1e-17, 1))
ylims!(axl, (1e-17, 1))
ylims!(axtl, (10^(-2.5), 10^(1.5)))
ylims!(axtr, (10^(-2.5), 10^(1.5)))

ks = range(0, 1.1 * π, length = 1000)

scrs = []
bs = []

for i in 1:3
    scatter!(axl, n_atoms, df_cube[i].error_rel, markersize = ms, marker = markers[i], color = colors[i])
    scr = scatter!(axr, n_atoms, df_thin[i].error_rel, markersize = ms, marker = markers[i], color = colors[i])
    
    hlines!(axl, [mean(df_cube[i].error_rel)], linestyle = :dash, linewidth = lw, color = colors[i])
    hlines!(axr, [mean(df_thin[i].error_rel)], linestyle = :dash, linewidth = lw, color = colors[i])
end

for i in 1:3
    preset = 2 * i
    mask = df_t_cube.preset .== preset
    b = round(FastSpecSoG.preset_parameters[2 * i][1], digits = 2)
    scatter!(axtl, n_atoms, df_t_cube.t_total[mask], markersize = ms, marker = markers[i], color = colors[i], label = L"%$b")

    @. model(x, p) = x + p[1]
    x_data = n_atoms
    y_data = df_t_cube.t_total[mask]
    p0 = [1.0]
    fit = curve_fit(model, log.(x_data), log.(y_data), p0)

    Ns = 10 .^ range(3, 5, 100)

    @info "(c) $i  x + $(fit.param[1])"
    g = x -> model(log(x), fit.param)
    lines!(axtl, Ns, exp.(g.(Ns)), linestyle = :dash, linewidth = lw, color = colors[i])
end

for i in 1:3
    preset = 2 * i
    mask = df_t_thin.preset .== preset
    scatter!(axtr, n_atoms, df_t_thin.t_total[mask], markersize = ms, marker = markers[i], color = colors[i])

    @. model(x, p) = x + p[1]
    x_data = n_atoms
    y_data = df_t_thin.t_total[mask]
    p0 = [1.0]
    fit = curve_fit(model, log.(x_data), log.(y_data), p0)

    Ns = 10 .^ range(3, 5, 100)

    @info "(d) $i x + $(fit.param[1])"
    g = x -> model(log(x), fit.param)
    lines!(axtr, Ns, exp.(g.(Ns)), linestyle = :dash, linewidth = lw, color = colors[i])
end

axislegend(axtl, L"b", position = :lt, nbanks = 1)

text!(axl, (1e5, 10^(-16.5)), text = "(a)", fontsize = 20, align = (:right, :baseline),)
text!(axr, (1e5, 10^(-16.5)), text = "(b)", fontsize = 20, align = (:right, :baseline),)
text!(axtl, (1e5, 10^(-2.4)), text = "(c)", fontsize = 20, align = (:right, :baseline),)
text!(axtr, (1e5, 10^(-2.4)), text = "(d)", fontsize = 20, align = (:right, :baseline),)

save("figs/total_error.pdf", f)
save("figs/total_error.png", f)

f