using ExTinyMD, FastSpecSoG, EwaldSummations, JLD2
using CSV, DataFrames

L = 10.0
N_real = (100, 100, 100)

w = (12, 12, 12)
β = 5.0 .* w
target = 2.0
#extra_pad_ratio_intial = ((target-1.0)*N_real0[1]+6.8*(target-1.0)*w[1]+2.4*target-4.4)/(2*w[1])
extra_pad_ratio = 12
lambda_true = (2*extra_pad_ratio*w[1]+N_real[1]+6.8*w[1]+4.4)/(N_real[1]+6.8*w[1]+2.4)
cheb_order = 16
preset = 6
uspara = USeriesPara(preset, r_c=1.799918355128706)
M_mid_initial = 18
eta = uspara.sw[M_mid_initial][1] / L + 0.0001
N_grid = (6, 6, 32)
Qs = 32
Ql = 32
r_c = 1.799918355128706
Q_0 = 32
Rz_0 = 32


path = "./reference/cube/n_1000_rc.jld2" 
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

for α in [1 / s0] ∪ [2.9:0.1:3.3...]
    s = r_c * α
    Ewald2D_interaction = Ewald2DInteraction(n_atoms, s, α, (L, L, L), ϵ = 1.0)
    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_e2d = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)
    err = abs(energy_ewald - energy_e2d)
    relerr = abs(energy_ewald - energy_e2d) / abs(energy_ewald)
    err_t = exp(-s^2) / s^2
    @show α, energy_e2d, err, relerr, err_t
end