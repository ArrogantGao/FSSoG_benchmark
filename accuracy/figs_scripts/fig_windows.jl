using CSV, DataFrames, CairoMakie, LaTeXStrings, SpecialFunctions, LsqFit

include("setting.jl")

df_ES = CSV.read("data/Acc_T3_1_WinFunc_ES.csv", DataFrame)
df_Gaussian = CSV.read("data/Acc_T3_1_WinFunc_Gaussian.csv", DataFrame)
df_PKB = CSV.read("data/Acc_T3_1_WinFunc_PKB.csv", DataFrame)

df_ESN = CSV.read("data/Acc_T3_2_Nspec_ES.csv", DataFrame)
df_GaussianN = CSV.read("data/Acc_T3_2_Nspec_Gaussian.csv", DataFrame)
df_PKBN = CSV.read("data/Acc_T3_2_Nspec_PKB.csv", DataFrame)

labels = ["KB", "ES", "Gaussian"]

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

axl = Axis(ga[1, 1], xlabel = L"\mathcal{P}", ylabel = L"\mathcal{E}_r", yscale = log10)
axr = Axis(gb[1, 1], xlabel = L"I", yscale = log10)

xlims!(axl, (0, 30))
xlims!(axr, (8, 68))
ylims!(axl, (1e-17, 1.0))
ylims!(axr, (1e-17, 1.0))

for (i, (l, error, m)) in enumerate(zip(labels, [error_PKB, error_ES, error_Gaussian], markers))
    scatter!(axl, w, error, markersize = ms, label = l, marker = m, color = colors[i])

    @. model(x, p) = log.(abs.(p[1] * erfc(p[2] * sqrt(x))))

    xs = [0.1:0.1:28.0...]

    raw_x_data = w
    raw_y_data = error
    mask_2 = abs.(raw_y_data) .> 1e-13
    x_data = raw_x_data[mask_2]
    y_data = raw_y_data[mask_2]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, x_data, log.(y_data), p0)

    @info "(a) $i $(fit.param[1]) erfc($(fit.param[2]) * sqrt(x))"

    g = x -> model(x, fit.param)
    lines!(axl, xs, exp.(g.(xs)), linestyle = :dash, linewidth = lw, color = colors[i])
end

for (i, (l, error, m)) in enumerate(zip(labels, [error_PKBN, error_ESN, error_GaussianN], markers))
    scatter!(axr, N, error, markersize = ms, label = l, marker = m, color = colors[i])

    #@. model(x, p) = p[1] * x + p[2]
    @. model(x, p) = log.(abs(p[1] * exp(- p[2]* x^2)))

    xs = [1:1:65...]

    raw_x_data = N
    raw_y_data = error
    mask_2 = abs.(raw_y_data) .> 1e-14
    x_data = raw_x_data[mask_2]
    y_data = raw_y_data[mask_2]
    p0 = [1.0, 0.01]
    fit = curve_fit(model, x_data, log.(y_data), p0)

    @info "(b) $i $(fit.param[1]) exp(-x^2 $(fit.param[2]))"

    g = x -> model(x, fit.param)
    lines!(axr, xs, exp.(g.(xs)), linestyle = :dash, linewidth = lw, color = colors[i])
end

axislegend(axl, position = :rt)

text!(axl, (29, 10^(-16.5)), text = "(a)", fontsize = 20, align = (:right, :baseline),)
text!(axr, (66, 10^(-16.5)), text = "(b)", fontsize = 20, align = (:right, :baseline),)

text!(axl, (30, 10^(-7.8)), text = L"O(\text{erfc}(C_1 \mathcal{P}^{-0.5}))", fontsize = 20, align = (:right, :baseline),)
text!(axr, (63, 10^(-7.8)), text = L"O(\exp(-C_2 I^{2}))", fontsize = 20, align = (:right, :baseline),)

save("figs/mid_windows.pdf", f)
save("figs/mid_windows.png", f)

f