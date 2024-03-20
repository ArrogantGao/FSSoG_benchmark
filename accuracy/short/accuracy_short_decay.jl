using FastSpecSoG, Plots, CSV, DataFrames

rs = [0.5:0.05:20.0...]

exacts = []

for preset in 1:5
    uspara = USeriesPara(preset)
    f_direct = x -> 1 / x - U_series(x, uspara)
    exact_values = f_direct.(rs)
    push!(exacts, exact_values)
end

plot(xlabel = "r", ylabel = "Error", title = "short range decay", dpi = 500)
for i in 1:5
    plot!(rs, abs.(exacts[i])
    , yscale = :log10, label = "preset = $i")
end
savefig("figs/accuracy_short_naive.png")