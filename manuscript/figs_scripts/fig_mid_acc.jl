using CSV, DataFrames, CairoMakie, LaTeXStrings, FastSpecSoG

df_1 = CSV.read("data/Acc_T2_1_Mmid.csv", DataFrame)
df_2 = CSV.read("data/Acc_T2_2_Extrapadding.csv", DataFrame)

uspara = USeriesPara(6)

λ_1 = unique(df_1.extra_pad_ratio)
λ_2 = unique(df_2.extra_pad_ratio)

Mmid_1 = df_1.M_mid
Mmid_2 = unique(df_2.M_mid)

marker = [:circle, :diamond, :star5, :utriangle]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axl = Axis(ga[1, 1], xlabel = L"s_{M}", ylabel = L"\mathcal{E}_r", yscale = log10, xscale = log10)
axr = Axis(gb[1, 1], xlabel = L"\lambda", yscale = log10)

ylims!(axl, (1e-16, 1))
ylims!(axr, (1e-16, 1))

sms = [uspara.sw[i][1] for i in 1:length(uspara.sw)]

for (i, λ) in enumerate(λ_1)
    mask = df_1.extra_pad_ratio .== λ
    scatter!(axl, sms[Mmid_1[mask]], df_1.error_rel[mask], markersize = 15, label = L"\lambda = %$λ", marker = marker[i])
end

for (i, Mmid) in enumerate(Mmid_2)
    mask = df_2.M_mid .== Mmid
    scatter!(axr, λ_2, df_2.error_rel[mask], markersize = 15, label = L"$M_{\text{mid}}$ = %$Mmid", marker = marker[i])
end

axislegend(axl, position = :lt)
axislegend(axr, position = :rt)

text!(axl, (100, 1e-13), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (30, 1e-13), text = "(b)", fontsize = 30, align = (:right, :baseline),)
f

save("figs/mid_acc.pdf", f)