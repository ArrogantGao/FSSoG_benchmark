using CSV, DataFrames, CairoMakie, LaTeXStrings, FastSpecSoG, SpecialFunctions, LsqFit

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
axr = Axis(gb[1, 1], xlabel = L"\lambda_z", yscale = log10)

ylims!(axl, (1e-16, 1))
ylims!(axr, (1e-16, 1))

sms = [uspara.sw[i][1] for i in 1:length(uspara.sw)]

pad_ratio = round.([138, 162, 186, 282] ./ (64 + 24), digits = 1)

for (i, λ) in enumerate(λ_1)
    mask = df_1.extra_pad_ratio .== λ
    λ_r = pad_ratio[i]
    scatter!(axl, sms[Mmid_1[mask]], df_1.error_rel[mask], markersize = 10, label = L"\lambda_z \approx %$(λ_r)", marker = marker[i])

    @. model(x, p) = log.(abs(p[1] * erfc(p[2] / x) ))

    xs = sms[Mmid_1[mask]]

    raw_x_data = sms[Mmid_1[mask]]
    raw_y_data = df_1.error_rel[mask]
    n0 = 1
    mask_2 = abs.(raw_y_data) .> 1e-14
    x_data = raw_x_data[mask_2]
    y_data = raw_y_data[mask_2]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, x_data, log.(y_data), p0)
    g = x -> model(x, fit.param)
    lines!(axl, xs, exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

f_pad = x -> (64 + 24 + 2 + x * 24) / (64 + 24)

for (i, Mmid) in enumerate(Mmid_2)
    mask = df_2.M_mid .== Mmid
    scatter!(axr, f_pad.(λ_2), df_2.error_rel[mask], markersize = 10, label = L"$M_{\text{mid}}$ = %$Mmid", marker = marker[i])

    @. model(x, p) = log.(abs(p[1] * erfc(p[2] * x) ))

    xs = λ_2

    raw_x_data = λ_2
    raw_y_data = df_2.error_rel[mask]
    mask_2 = abs.(raw_y_data) .> 1e-14
    x_data = raw_x_data[mask_2]
    y_data = raw_y_data[mask_2]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, x_data, log.(y_data), p0)
    g = x -> model(x, fit.param)
    lines!(axr, f_pad.(xs), exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

axislegend(axl, position = :lt)
axislegend(axr, position = :rt)

text!(axl, (100, 1e-13), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (9.3, 1e-13), text = "(b)", fontsize = 30, align = (:right, :baseline),)

text!(axl, (100, 10^(-9.5)), text = L"O(\text{erfc}(C_1 s_M^{-1}))", fontsize = 20, align = (:right, :baseline),)
text!(axr, (9.6, 10^(-9.5)), text = L"O(\text{erfc}(C_2 \lambda_z))", fontsize = 20, align = (:right, :baseline),)

save("figs/mid_acc.pdf", f)
save("figs/mid_acc.png", f)

f