using FastSpecSoG, Plots, CSV, DataFrames, Random, ExTinyMD, LaTeXStrings
Random.seed!(1234)

n_atoms0 = 1000
N_grid0 = (40, 40, 40)
N_grid_exact0 = (80, 80, 80)
L0 = (20.0, 20.0, 20.0)
ratios = [2.0 ^i for i in 0.0:0.5:3.0]

w = (8, 8, 8)
β = 5.0 .* w
extra_pad_ratios = [2, 3, 4, 8]
cheb_order = 10

w_exact = (16, 16, 16)
β_exact = 5.0 .* w_exact
extra_pad_ratio_exact = 8
cheb_order_exact = 24

uspara = USeriesPara(6)
M_mid_initial = 14
eta = uspara.sw[M_mid_initial][1] / L0[3]


errors = [[] for i in 1:length(extra_pad_ratios)]
Num_atoms = [Int(ceil(n_atoms0 * ratio^3)) for ratio in ratios]

for ratio in ratios
    n_atoms = Int(ceil(n_atoms0 * ratio^3))
    N_grid = Int.(ceil.(N_grid0 .* ratio))
    N_grid_exact = Int.(ceil.(N_grid_exact0 .* ratio))
    L = L0 .* ratio

    M_mid = M_mid_initial
    
    #for i in 1:length(uspara.sw)
        #if (uspara.sw[i][1] < eta * L[3]) && (uspara.sw[i+1][1] > eta * L[3])
            #M_mid = i
            #break
        #end
    #end

    @show n_atoms, N_grid, N_grid_exact, L, M_mid

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
        CSV.write("data/Acc_T3_1_window.csv", df, append = true)
    end
end
