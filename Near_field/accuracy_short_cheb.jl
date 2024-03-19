using FastSpecSoG, Plots

terms = [1:3:40...]
errors = [[] for i in 1:5]

r_min = 0.5
r_max = 10.0
rs = [r_min:0.01:r_max...]

for preset in 1:5
    uspara = USeriesPara(preset)
    f_direct = x -> U_series(x, uspara)
    exact_values = f_direct.(rs)
    for n in terms
        uspara_cheb = FastSpecSoG.Es_USeries_Cheb(uspara, r_min, r_max, n)
        f_cheb = x -> uspara_cheb(x)
        cheb_values = f_cheb.(rs)
        push!(errors[preset], maximum(abs.(cheb_values .- exact_values)))
    end
end

plot(xlabel = "Cheb order", ylabel = "Max error", title = "Accuracy of Chebyshev series", dpi = 500)
for i in 1:5
    plot!(terms, errors[i], label = "Preset $i", yscale = :log10, yaxis = (1e-16, 1.0), marker = :circle)
end
display(plot!())

savefig("Near_field/figs/accuracy_short_cheb.png")