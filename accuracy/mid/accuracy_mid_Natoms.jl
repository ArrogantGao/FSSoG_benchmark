using FastSpecSoG, Plots, CSV, DataFrames, Random, ExTinyMD, LaTeXStrings
Random.seed!(1234)

n_atoms0 = 100
N_grid0 = (64, 64, 64)
N_grid_exact0 = (128, 128, 128)
L0 = (100.0, 100.0, 100.0)
ratios = [2.0 ^i for i in 0.0:0.25:2.0]

w = (16, 16, 16)
β = 5.0 .* w
extra_pad_ratios = [2, 3, 4, 8]
cheb_order = 10

w_exact = (32, 32, 32)
β_exact = 5.0 .* w_exact
extra_pad_ratio_exact = 8
cheb_order_exact = 24

uspara = USeriesPara(2)
M_mid = 6

errors = [[] for i in 1:length(extra_pad_ratios)]
Num_atoms = [Int(ceil(n_atoms0 * ratio^3)) for ratio in ratios]

for ratio in ratios
    n_atoms = Int(ceil(n_atoms0 * ratio^3))
    N_grid = Int.(ceil.(N_grid0 .* ratio))
    N_grid_exact = Int.(ceil.(N_grid_exact0 .* ratio))
    L = L0 .* ratio

    @show n_atoms, N_grid, N_grid_exact, L

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_grid_exact, w_exact, β_exact, L, extra_pad_ratio_exact, cheb_order_exact, uspara, M_mid)
    @time E_exact = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

    for extra_pad_ratio in extra_pad_ratios

        gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_grid, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
        @time E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

        abs_error = abs(E_3DFFT - E_exact)
        relative_error = abs_error / abs(E_exact)
        @show extra_pad_ratio, abs_error, relative_error

        df = DataFrame(extra_pad_ratio = extra_pad_ratio, N_atoms = n_atoms, energy_exact = E_exact, energy_approx = E_3DFFT, abs_error = abs_error, relative_error = relative_error)
        CSV.write("data/accuracy_mid_Natoms.csv", df, append = true)
    end
end

plot(xlabel = L"N_{atoms}", ylabel = L"Error", title = L"Error vs Number of atoms", dpi = 500, xlim = [2, 4])
for i in 1:length(extra_pad_ratios)
    plot!(log10.(Num_atoms), log10.(abs.(errors[i])), label = "λ = $(extra_pad_ratios[i])", marker = :circle)
end
savefig("figs/accuracy_mid_Natoms.png")