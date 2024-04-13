using CSV, DataFrames, CairoMakie, LaTeXStrings, FastSpecSoG, LsqFit

df_cube_2 = CSV.read("data/Acc_T6_cube_total_2.csv", DataFrame)
df_cube_4 = CSV.read("data/Acc_T6_cube_total_4.csv", DataFrame)
df_cube_6 = CSV.read("data/Acc_T6_cube_total_6.csv", DataFrame)
df_thin_2 = CSV.read("data/Acc_T6_thin_total_2.csv", DataFrame)
df_thin_4 = CSV.read("data/Acc_T6_thin_total_4.csv", DataFrame)
df_thin_6 = CSV.read("data/Acc_T6_thin_total_6.csv", DataFrame)

df_cube = [df_cube_2, df_cube_4, df_cube_6]
df_thin = [df_thin_2, df_thin_4, df_thin_6]

n_atoms = [1000, 3164, 10000, 31624, 100000]

marker = [:circle, :diamond, :star5, :utriangle, :hexagon, :xcross]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axr = Axis(gb[1, 1], xlabel = L"N_\text{atoms}", yscale = log10, xscale = log10)
axl = Axis(ga[1, 1], xlabel = L"N_\text{atoms}", ylabel = L"\mathcal{E}_r", yscale = log10, xscale = log10)

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

axislegend(axr, scrs, bs, position = :rb, nbanks = 2)

text!(axl, (2000, 1e-2), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (2000, 1e-2), text = "(b)", fontsize = 30, align = (:right, :baseline),)

save("figs/total_error.pdf", f)
save("figs/total_error.png", f)

f