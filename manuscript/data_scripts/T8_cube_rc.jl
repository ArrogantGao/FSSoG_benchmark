using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

L = 10.0
N_real = (94, 94, 94)

w = (7, 7, 7)
β = 5.0 .* w
target = 1.5
#extra_pad_ratio_intial = ((target-1.0)*N_real0[1]+6.8*(target-1.0)*w[1]+2.4*target-4.4)/(2*w[1])
extra_pad_ratio = 10
lambda_true = (2*extra_pad_ratio*w[1]+N_real[1]+6.8*w[1]+4.4)/(N_real[1]+6.8*w[1]+2.4)
@show lambda_true
cheb_order = 10
preset = 6
uspara = USeriesPara(preset, r_c=1.799918355128706)
M_mid_initial = 18
eta = uspara.sw[M_mid_initial][1] / L + 0.0001
@show eta
N_grid = (3, 3, 20)
Qs = 25
Ql = 18
r_c = 1.799918355128706
Q_0 = 20
Rz_0 = 20


path = "./reference/cube/n_100_rc.jld2" 
@load path n_atoms L atoms info energy_ewald 

@show n_atoms, L, energy_ewald

boundary = ExTinyMD.Q2dBoundary(L, L, L)

M_mid = proper_M(eta, L, uspara)
fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Qs, 0.1, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Ql, Rz_0, Q_0; preset = preset, ϵ = 1.0)
fssog_neighbor = CellList3D(info, fssog_interaction.r_c, boundary, 1)
@time energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

abs_error = abs(energy_ewald - energy_sog)
relative_error = abs(energy_ewald - energy_sog) / abs(energy_ewald)

@show energy_sog, abs_error, relative_error
s0 = fssog_interaction.uspara.sw[1][1]
@show 1 / s0

for α in 0.1:0.1:3.0
    s = r_c * α
    Ewald2D_interaction = Ewald2DInteraction(n_atoms, s, α, (L, L, L), ϵ = 1.0)
    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_e2d = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)
    err = abs(energy_ewald - energy_e2d)
    relerr = abs(energy_ewald - energy_e2d) / abs(energy_ewald)
    @show α, energy_e2d, err, relerr
end