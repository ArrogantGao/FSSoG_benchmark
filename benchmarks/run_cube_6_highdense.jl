using FastSpecSoG, ExTinyMD, BenchmarkTools, Plots, CSV, DataFrames, ArgParse, JLD2

function parse_commandline()

  s = ArgParseSettings()

  @add_arg_table s begin
    "--output-file"
    help = "relative file path of the output file"
  end

  return parse_args(s)

end

function main()
    
    parsed_args = parse_commandline()
    output_file = parsed_args["output-file"]


    L0 = 10.0
N_real0 = (94, 94, 94)

w = (7, 7, 7)
β = 5.0 .* w
target = 1.5
#extra_pad_ratio_intial = ((target-1.0)*N_real0[1]+6.8*(target-1.0)*w[1]+2.4*target-4.4)/(2*w[1])
extra_pad_ratio_intial = 10
lambda_true = (2*extra_pad_ratio_intial*w[1]+N_real0[1]+6.8*w[1]+4.4)/(N_real0[1]+6.8*w[1]+2.4)
@show lambda_true
cheb_order = 10
preset = 6
uspara = USeriesPara(preset, r_c=1.799918355128706)
M_mid_initial = 18
eta = uspara.sw[M_mid_initial][1] / L0 + 0.0001
@show eta
N_grid = (3, 3, 20)
Qs = 25
Ql = 18
r_c = 1.799918355128706
Q_0 = 20
Rz_0 = 20

    for data in ["n_100000.jld2"]

        path = "../manuscript/reference/cube/high_rho/$data" 
        @load path n_atoms L atoms info energy_ewald

        boundary = ExTinyMD.Q2dBoundary(L, L, L)

        ratio = (n_atoms / 100000)^(1/3)
        N_real = Int.(ceil.(ratio .* N_real0))
        extra_pad_ratio = Int(ceil(extra_pad_ratio_intial.*ratio))

        M_mid = proper_M(eta, L, uspara)
        interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Qs, 0.1, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Ql, Rz_0, Q_0; preset = preset, ϵ = 1.0)
        neighbor = CellList3D(info, interaction.r_c, boundary, 1)
        

        for i in 1:interaction.n_atoms
            interaction.charge[i] = atoms[info.particle_info[i].id].charge
            interaction.position[i] = info.particle_info[i].position.coo
        end

        t_short = @belapsed energy_short($interaction, $neighbor)
        t_mid = @belapsed energy_mid($interaction)
        t_long = @belapsed energy_long($interaction)

        @show preset, n_atoms, t_short, t_mid, t_long

        df = DataFrame(preset = preset, n_atoms = n_atoms, t_short = t_short, t_mid = t_mid, t_long = t_long, t_total = t_short + t_mid + t_long)
        CSV.write(output_file, df, append = true)
    end
end

main()