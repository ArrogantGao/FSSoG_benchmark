using CSV, DataFrames, CairoMakie, LaTeXStrings, FastSpecSoG, LsqFit

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

marker = [:circle, :diamond, :star5, :utriangle, :hexagon, :xcross]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 700), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

gta = f[2, 1] = GridLayout()
gtb = f[2, 2] = GridLayout()

axr = Axis(gb[1, 1], xlabel = L"N_\text{atoms}", yscale = log10, xscale = log10)
axl = Axis(ga[1, 1], xlabel = L"N_\text{atoms}", ylabel = L"\mathcal{E}_r", yscale = log10, xscale = log10)

axtl = Axis(gta[1, 1], xlabel = L"N_\text{atoms}", ylabel = L"\text{Time (s)}", yscale = log10, xscale = log10)
axtr = Axis(gtb[1, 1], xlabel = L"N_\text{atoms}", yscale = log10, xscale = log10)

ylims!(axr, (1e-16, 1))
ylims!(axl, (1e-16, 1))

ks = range(0, 1.1 * π, length = 1000)

scrs = []
bs = []

for i in 1:3
    b = round(FastSpecSoG.preset_parameters[2 * i][1], digits = 2)
    scatter!(axl, n_atoms, df_cube[i].error_rel, markersize = 10, marker = marker[i])
    scr = scatter!(axr, n_atoms, df_thin[i].error_rel, markersize = 10, marker = marker[i])
    push!(bs, L"b = %$b")
    push!(scrs, scr)
end

for i in 1:3
    preset = 2 * i
    mask = df_t_cube.preset .== preset
    scatter!(axtl, n_atoms, df_t_cube.t_total[mask], markersize = 10, marker = marker[i])

    @. model(x, p) = x + p[1]
    x_data = n_atoms
    y_data = df_t_cube.t_total[mask]
    p0 = [1.0]
    fit = curve_fit(model, log.(x_data), log.(y_data), p0)

    Ns = 10 .^ range(3, 5, 100)

    @info "(c) $i  x + $(fit.param[1])"
    g = x -> model(log(x), fit.param)
    lines!(axtl, Ns, exp.(g.(Ns)), linestyle = :dash, linewidth = 0.7)
end

for i in 1:3
    preset = 2 * i
    mask = df_t_thin.preset .== preset
    scatter!(axtr, n_atoms, df_t_thin.t_total[mask], markersize = 10, marker = marker[i])

    @. model(x, p) = x + p[1]
    x_data = n_atoms
    y_data = df_t_thin.t_total[mask]
    p0 = [1.0]
    fit = curve_fit(model, log.(x_data), log.(y_data), p0)

    Ns = 10 .^ range(3, 5, 100)

    @info "(d) $i x + $(fit.param[1])"
    g = x -> model(log(x), fit.param)
    lines!(axtr, Ns, exp.(g.(Ns)), linestyle = :dash, linewidth = 0.7)
end

axislegend(axtl, scrs, bs, position = :rb, nbanks = 1)

text!(axl, (2000, 1e-2), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (2000, 1e-2), text = "(b)", fontsize = 30, align = (:right, :baseline),)
text!(axtl, (2000, 5), text = "(c)", fontsize = 30, align = (:right, :baseline),)
text!(axtr, (2000, 10), text = "(d)", fontsize = 30, align = (:right, :baseline),)

save("figs/total_error.pdf", f)
save("figs/total_error.png", f)

f