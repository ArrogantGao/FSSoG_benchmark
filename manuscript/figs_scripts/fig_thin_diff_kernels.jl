using CSV, DataFrames, CairoMakie, LaTeXStrings, LsqFit

df_w = CSV.read("data/Acc_T5_thin_w.csv", DataFrame)
df_Taylor = CSV.read("data/Acc_T5_thin_TaylorQ.csv", DataFrame)
df_Nxy = CSV.read("data/Acc_T5_thin_Nxy.csv", DataFrame)
df_Rz = CSV.read("data/Acc_T5_thin_Rz.csv", DataFrame)

markers = [:circle, :diamond, :star5, :utriangle]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 800), fontsize = 20)
ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()
gc = f[2, 1] = GridLayout()
gd = f[2, 2] = GridLayout()

ax_w = Axis(ga[1, 1], xlabel = L"w", ylabel = L"\mathcal{E}_r", yscale = log10)
ax_Taylor = Axis(gb[1, 1], xlabel = L"Q_{\text{Taylor}}", yscale = log10)
ax_Nxy = Axis(gc[1, 1], xlabel = L"N_{xy}", ylabel = L"\mathcal{E}_r", yscale = log10, xscale = log2)
ax_Rz = Axis(gd[1, 1], xlabel = L"R_z", yscale = log10)

xlims!(ax_w, (0, 17))
xlims!(ax_Taylor, (0, 17))
xlims!(ax_Nxy, (2^3.5, 2^10.5))
xlims!(ax_Rz, (0, 17))

ylims!(ax_w, (1e-16, 1))
ylims!(ax_Taylor, (1e-16, 1))
ylims!(ax_Nxy, (1e-16, 1))
ylims!(ax_Rz, (1e-16, 1))

df_w_errors = [df_w.error_rel_pkb, df_w.error_ES, df_w.error_Gaussian]
labels = ["PKB", "ES", "Gaussian"]
for i in 1:3
    scatter!(ax_w, df_w.w, df_w_errors[i], markersize = 10, label = labels[i], marker = markers[i])
end

for i in 1:1
    scatter!(ax_Taylor, df_Taylor.TaylorQ, df_Taylor.error_rel, markersize = 10, marker = markers[i])
    scatter!(ax_Nxy, df_Nxy.Nxy, df_Nxy.error_rel, markersize = 10, marker = markers[i])
    scatter!(ax_Rz, df_Rz.R_z, df_Rz.error_rel, markersize = 10, marker = markers[i])
end

text!(ax_w, (4, 1e-13), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Taylor, (4, 1e-13), text = "(b)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Nxy, (2^5, 1e-13), text = "(c)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Rz, (4, 1e-13), text = "(d)", fontsize = 30, align = (:right, :baseline),)

axislegend(ax_w, position = :rt)

save("figs/fig_thin.pdf", f)
save("figs/fig_thin.png", f)

f