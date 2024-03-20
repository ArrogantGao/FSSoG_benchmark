using FastSpecSoG, StatProfilerHTML, BenchmarkTools

n_atoms = 1000
L = (100.0, 100.0, 100.0)

qs = rand(n_atoms)
qs .-= sum(qs) ./ n_atoms
poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

# using the 3DFFT
N_real = (16, 16, 16)
w = Int.(N_real ./ 4)
β = 5.0 .* w
extra_pad_ratio = 2
cheb_order = 10
uspara = USeriesPara(2)
M_mid = 3

gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)

energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

@show @belapsed energy_mid($qs, $poses, $gridinfo, $gridbox, $cheb_coefs, $scalefactor)

@profilehtml for _=1:1000
    energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)
end