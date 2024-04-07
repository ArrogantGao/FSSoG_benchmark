using FastSpecSoG, CSV, DataFrames
using Random
Random.seed!(123)

n_atoms_0 = 1000
L_0 = (20.0, 20.0, 20.0)
N_real_0 = (256, 256)

ratios = [10^i for i in 0:0.5:3]

for ratio in ratios
    n_atoms = Int(ceil(n_atoms_0 * ratio / 2) * 2)
    L = L_0 .* ratio^(1/3)
    qs = [(-1.0)^i for i in 1:n_atoms]
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    # exact results
    uspara = USeriesPara(6)
    M_total = length(uspara.sw)
    N_real = (128, 128, 128)
    


end