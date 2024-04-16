using CSV, DataFrames, CairoMakie, LaTeXStrings, LsqFit, SpecialFunctions

df_w = CSV.read("data/Acc_T5_thin_w.csv", DataFrame)
df_Taylor = CSV.read("data/Acc_T5_thin_TaylorQ.csv", DataFrame)
df_Nxy = CSV.read("data/Acc_T5_thin_Nxy.csv", DataFrame)
df_Rz = CSV.read("data/Acc_T5_thin_Rz.csv", DataFrame)

gammas = unique(df_w.gamma)

markers = [:circle, :diamond, :star5, :utriangle]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 700), fontsize = 20)
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

    @. model(x, p) = log.(abs.(p[1] * erfc(p[2] * sqrt(x))))

    xs = [0.1:0.1:16.0...]

    raw_x_data = df_w.w[mask_w]
    raw_y_data = df_w.error_rel[mask_w]
    mask_2 = abs.(raw_y_data) .> 1e-13
    x_data = raw_x_data[mask_2]
    y_data = raw_y_data[mask_2]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, x_data, log.(y_data), p0)

    @info "(a) $i $(fit.param[1]) erfc($(fit.param[2]) * sqrt(x))"

    g = x -> model(x, fit.param)
    lines!(ax_w, xs, exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

for i in 1:3
    gamma = gammas[i]
    mask_Taylor = df_Taylor.gamma .== gamma
    scatter!(ax_Taylor, df_Taylor.TaylorQ[mask_Taylor], df_Taylor.error_rel[mask_Taylor], markersize = 10, marker = markers[i])

    @. model(x, p) = log(abs(p[2] / (factorial(x) * (p[1])^(2x) )))

    xs = [0:16...]

    raw_x_data = df_Taylor.TaylorQ[mask_Taylor]
    raw_y_data = df_Taylor.error_rel[mask_Taylor]
    mask_2 = abs.(raw_y_data) .> 1e-12
    x_data = raw_x_data[mask_2]
    y_data = raw_y_data[mask_2]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, (x_data), log.(y_data), p0)

    @info "(b) $i $(fit.param[2]) / (x! * (fit.param[1])^(2x) ))"

    g = x -> model(x, fit.param)
    lines!(ax_Taylor, xs, exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

for i in 1:3
    gamma = gammas[i]
    mask_Nxy = df_Nxy.gamma .== gamma
    scatter!(ax_Nxy, df_Nxy.Nxy[mask_Nxy], df_Nxy.error_rel[mask_Nxy], markersize = 10, marker = markers[i])

    @. model(x, p) = p[1] * x + p[2]

    xs = range(2^3.5, 2^10.5, length = 1000)

    raw_x_data = df_Nxy.Nxy[mask_Nxy]
    raw_y_data = df_Nxy.error_rel[mask_Nxy]
    mask_2 = abs.(raw_y_data) .> 1e-12
    x_data = raw_x_data[mask_2]
    y_data = raw_y_data[mask_2]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, (x_data), log.(y_data), p0)

    @info "(c) $i $(fit.param[1]) x + $(fit.param[2])"

    g = x -> model(x, fit.param)
    lines!(ax_Nxy, xs, exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

for i in 1:3
    gamma = gammas[i]
    mask_Rz = df_Rz.gamma .== gamma
    scatter!(ax_Rz, df_Rz.R_z[mask_Rz], df_Rz.error_rel[mask_Rz], markersize = 10, marker = markers[i])

    @. model(x, p) = log.(abs(p[1] / (p[2]^x * sqrt(factorial(big(x)))) ))

    xs = [0:15...]

    raw_x_data = df_Rz.R_z[mask_Rz]
    raw_y_data = df_Rz.error_rel[mask_Rz]
    mask_2 = abs.(raw_y_data) .> 1e-13
    x_data = raw_x_data[mask_2]
    y_data = raw_y_data[mask_2]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, x_data, log.(y_data), p0)

    @info "(d) $i $(fit.param[1]) / ( $(fit.param[2])^x * sqrt(x!))"

    g = x -> model(x, fit.param)
    lines!(ax_Rz, xs, exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

text!(ax_w, (15.5, 10^(-7.8)), text = L"O(\text{erfc}(C_1 w^{-0.5}))", fontsize = 20, align = (:right, :baseline),)
text!(ax_Nxy, (2^10, 10^(-5.8)), text = L"O(\exp( - C_2 N_{xy}))", fontsize = 20, align = (:right, :baseline),)
text!(ax_Taylor, (12, 10^(-7.8)), text = L"O((Q! C_3^{2Q})^{-1})", fontsize = 20, align = (:right, :baseline),)
text!(ax_Rz, (16, 10^(-7.8)), text = L"O(C_4^{-R_z} / \sqrt{R_z !})", fontsize = 20, align = (:right, :baseline),)

text!(ax_w, (2.5, 1e-13), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Nxy, (2^4.5, 1e-13), text = "(b)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Taylor, (2.5, 1e-13), text = "(c)", fontsize = 30, align = (:right, :baseline),)
text!(ax_Rz, (2.5, 1e-13), text = "(d)", fontsize = 30, align = (:right, :baseline),)

axislegend(ax_w, position = :rt)

save("figs/fig_thin.pdf", f)
save("figs/fig_thin.png", f)

f