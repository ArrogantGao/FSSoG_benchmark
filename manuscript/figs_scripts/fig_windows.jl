using CSV, DataFrames, CairoMakie, LaTeXStrings

df_ES = CSV.read("data/Acc_T3_1_WinFunc_ES.csv", DataFrame)
df_Gaussian = CSV.read("data/Acc_T3_1_WinFunc_Gaussian.csv", DataFrame)
df_PKB = CSV.read("data/Acc_T3_1_WinFunc_PKB.csv", DataFrame)

df_ESN = CSV.read("data/Acc_T3_2_Nspec_ES.csv", DataFrame)
df_GaussianN = CSV.read("data/Acc_T3_2_Nspec_Gaussian.csv", DataFrame)
df_PKBN = CSV.read("data/Acc_T3_2_Nspec_PKB.csv", DataFrame)

labels = ["ES", "Gaussian", "PKB"]
marker = [:circle, :diamond, :star5]

w = df_ES.w
error_ES = df_ES.error_rel
error_Gaussian = df_Gaussian.error_rel
error_PKB = df_PKB.error_rel

N = df_ESN.N_real
error_ESN = df_ESN.error_rel
error_GaussianN = df_GaussianN.error_rel
error_PKBN = df_PKBN.error_rel

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axl = Axis(ga[1, 1], xlabel = L"w", ylabel = L"\mathcal{E}_r", yscale = log10)
axr = Axis(gb[1, 1], xlabel = L"N_{\text{real}}", yscale = log10)

xlims!(axl, (0, 13))
xlims!(axr, (8, 68))
ylims!(axl, (1e-16, 1.0))
ylims!(axr, (1e-16, 1.0))

for (l, error, m) in zip(labels, [error_ES, error_Gaussian, error_PKB], marker)
    scatter!(axl, w, error, markersize = 10, label = l, marker = m)
end

for (l, error, m) in zip(labels, [error_ESN, error_GaussianN, error_PKBN], marker)
    scatter!(axr, N, error, markersize = 10, label = l, marker = m)
end

axislegend(axl, position = :rt)

text!(axl, (3, 1e-14), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (22.5, 1e-14), text = "(b)", fontsize = 30, align = (:right, :baseline),)

save("figs/mid_windows.pdf", f)
save("figs/mid_windows.png", f)