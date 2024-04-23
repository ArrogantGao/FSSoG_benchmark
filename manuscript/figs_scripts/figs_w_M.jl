using FastSpecSoG, CSV, DataFrames
using CairoMakie, LaTeXStrings

N = 94

df_w_M = CSV.read("data/cube_w_M.csv", DataFrame)

M_mids = unique(df_w_M.M_mid)

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (400, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
axl = Axis(ga[1, 1], xlabel = L"w", ylabel = L"\mathcal{E}_r", yscale = log10, title = "N = ($N, $N, $N)")

ylims!(axl, 1e-16, 1)

for M_mid in M_mids
    mask = df_w_M.M_mid .== M_mid
    scatter!(axl, df_w_M.w[mask], abs.(df_w_M.rel_error[mask] .+ 1e-16), markersize = 10, label = L"M_1 = %$M_mid")
end

axislegend(axl, position = :rt)
save("figs/figs_w_M_$(N).pdf", f)
f