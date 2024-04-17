using CSV, DataFrames, CairoMakie, LaTeXStrings, LsqFit, FastSpecSoG

df_energy = CSV.read("data/Acc_T0_energy.csv", DataFrame)
df_force = CSV.read("data/Acc_T0_force.csv", DataFrame)

Ms = unique(df_force.M)
Mse = unique(df_energy.M)
presets = unique(df_force.preset)

marker = [:circle, :diamond, :star5, :utriangle, :hexagon, :xcross]

f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (800, 400), fontsize = 20)

ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axl = Axis(ga[1, 1], xlabel = L"M", ylabel = L"\mathcal{E}_r", yscale = log10)
axr = Axis(gb[1, 1], xlabel = L"M", yscale = log10)

xlims!(axl, (0, 300))
xlims!(axr, (0, 120))

ylims!(axl, (1e-16, 1))
ylims!(axr, (1e-16, 1))

global ss1 = 0.0
for preset in presets
    b0 = FastSpecSoG.preset_parameters[preset][1]
    c = log(b0)^(-1.5) * exp(-π^2 / (2 * log(b0)))
    mask_energy = df_energy.preset .== preset
    cr1 = df_energy.error[mask_energy][end]
    ss1 += cr1 / c
end
global as1 = ss1 / 6
@show as1


for preset in presets
    mask_energy = df_energy.preset .== preset
    b0 = FastSpecSoG.preset_parameters[preset][1]
    b = round(FastSpecSoG.preset_parameters[preset][1], digits = 2)
    scatter!(axl, Mse, df_energy.error[mask_energy], markersize = 10, label = L"b = %$b", marker = marker[preset])

    c = log(b0)^(-1.5) * exp(-π^2 / (2 * log(b0)))

    @. model(x, p) = log((c * as1 + 1 / b0^(x + p[1])))
    raw_y_data = df_energy.error[mask_energy]
    raw_x_data = Mse

    x_data = raw_x_data[1:end]
    y_data = raw_y_data[1:end]
    p0 = [1.0, 1.0]
    fit = curve_fit(model, x_data, log.(abs.(y_data)), p0)

    # @info "(a) $preset $(df_energy.error[mask_energy][end]) + $(fit.param[1]) / $b0^x"
    g = x -> model(x, fit.param)
    xs = [0.0:0.1:300.0...]
    lines!(axl, xs, exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

global ss = 0.0
for preset in presets
    b0 = FastSpecSoG.preset_parameters[preset][1]
    c = log(b0)^(-1.5) * exp(-π^2 / (2 * log(b0)))
    mask_force = df_force.preset .== preset
    cr = df_force.error[mask_force][end]
    ss += cr / c
end
global as = ss / 6
@show as

for preset in presets
    mask_force = df_force.preset .== preset
    b0 = FastSpecSoG.preset_parameters[preset][1]
    b = round(FastSpecSoG.preset_parameters[preset][1], digits = 2)
    scatter!(axr, Ms, df_force.error[mask_force], markersize = 10, label = L"b = %$b", marker = marker[preset])

    c = log(b0)^(-1.5) * exp(-π^2 / (2 * log(b0)))

    @. model(x, p) = log((c * as + 1 / b0^(3 * x + p[1])))
    raw_y_data = df_force.error[mask_force]
    raw_x_data = Ms

    x_data = raw_x_data[preset:end]
    y_data = raw_y_data[preset:end]
    p0 = [1.0]
    fit = curve_fit(model, x_data, log.(abs.(y_data)), p0)

    @show preset b0^fit.param[1]

    # @info "(b) $preset $(df_force.error[mask_force][end]) + $(fit.param[1]) / $b0^3x"
    g = x -> model(x, fit.param)
    xs = [0.0:0.1:120.0...]
    lines!(axr, xs, exp.(g.(xs)), linestyle = :dash, linewidth = 0.7)
end

axislegend(axl, position = :lb, nbanks = 2)

text!(axl, (50, 1e-11), text = "(a)", fontsize = 30, align = (:right, :baseline),)
text!(axr, (20, 1e-11), text = "(b)", fontsize = 30, align = (:right, :baseline),)

# text!(axl, (180, 1e-11), text = L"O(b^{-M} + C_1)", fontsize = 20, align = (:right, :baseline),)
# text!(axr, (72, 1e-11), text = L"O(b^{-3M} + C_2)", fontsize = 20, align = (:right, :baseline),)

f

save("figs/cutoff.pdf", f)
save("figs/cutoff.png", f)