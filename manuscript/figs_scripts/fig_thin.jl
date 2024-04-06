using CSV, DataFrames, CairoMakie, LaTeXStrings

df_w = CSV.read("data/Acc_T5_thin_w.csv", DataFrame)
df_Taylor = CSV.read("data/Acc_T5_thin_TaylorQ.csv", DataFrame)
df_Nxy = CSV.read("data/Acc_T5_thin_Nxy.csv", DataFrame)
df_Rz = CSV.read("data/Acc_T5_thin_Rz.csv", DataFrame)

gammas = unique(df_w.gamma)

markers = [:circle, :diamond, :star5, :utriangle]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 800), fontsize = 20)
ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()
gc = f[2, 1] = GridLayout()
gd = f[2, 2] = GridLayout()

ax_w = Axis(ga[1, 1], xlabel = L"w", ylabel = L"\mathcal{E}_r", yscale = log10)
ax_Nxy = Axis(gb[1, 1], xlabel = L"N_{xy}", yscale = log10, xscale = log2)
ax_Taylor = Axis(gc[1, 1], xlabel = L"Q_{\text{Taylor}}", ylabel = L"\mathcal{E}_r", yscale = log10)
ax_Rz = Axis(gd[1, 1], xlabel = L"R_z", yscale = log10)

xlims!(ax_w, (0, 17))
xlims!(ax_Taylor, (0, 17))
xlims!(ax_Nxy, (2^3.5, 2^10.5))
xlims!(ax_Rz, (0, 17))

ylims!(ax_w, (1e-16, 1))
ylims!(ax_Taylor, (1e-16, 1))
ylims!(ax_Nxy, (1e-16, 1))
ylims!(ax_Rz, (1e-16, 1))


for i in 1:3
    gamma = gammas[i]
    mask_w = df_w.gamma .== gamma
    scatter!(ax_w, df_w.w[mask_w], df_w.error_rel[mask_w], markersize = 10, label = L"\gamma = %$gamma", marker = markers[i])
end

for i in 1:3
    gamma = gammas[i]
    mask_Taylor = df_Taylor.gamma .== gamma
    scatter!(ax_Taylor, df_Taylor.TaylorQ[mask_Taylor], df_Taylor.error_rel[mask_Taylor], markersize = 10, marker = markers[i])
end

for i in 1:3
    gamma = gammas[i]
    mask_Nxy = df_Nxy.gamma .== gamma
    scatter!(ax_Nxy, df_Nxy.Nxy[mask_Nxy], df_Nxy.error_rel[mask_Nxy], markersize = 10, marker = markers[i])
end

for i in 1:3
    gamma = gammas[i]
    mask_Rz = df_Rz.gamma .== gamma
    scatter!(ax_Rz, df_Rz.R_z[mask_Rz], df_Rz.error_rel[mask_Rz], markersize = 10, marker = markers[i])
end

text!(ax_w, (2.5, 1e-13), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Nxy, (2^4.7, 1e-13), text = "(b)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Taylor, (2.5, 1e-13), text = "(c)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Rz, (2.5, 1e-13), text = "(d)", fontsize = 30, align = (:right, :baseline),)

axislegend(ax_w, position = :rt)

save("figs/fig_thin.pdf", f)
save("figs/fig_thin.png", f)

f