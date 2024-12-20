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


    L0 = 20.0

    N_real0 = (36, 36, 36)
    w = (8, 8, 8)
    β = 5.0 .* w
    extra_pad_ratio_intial = 7
    lambda_true = (2*extra_pad_ratio_intial*w[1]+N_real0[1]+2*w[1]+2)/(N_real0[1]+2*w[1])
    @show lambda_true
    cheb_order = 10
    preset = 6
    uspara = USeriesPara(preset) 
    M_mid_initial = 15
    eta = uspara.sw[M_mid_initial][1] / L0 + 0.0001
    @show eta
    N_grid = (2, 2, 14)
    Qs = 32
    Ql = 16
    r_c = 9.99
    Rz_0 = 15
    Q_0 = 20

    for data in ["n_1000.jld2", "n_3164.jld2", "n_10000.jld2", "n_31624.jld2", "n_100000.jld2", "n_316228.jld2", "n_1000000.jld2"]

      path = "../manuscript/reference/cube_200/$data"
      @load path n_atoms L atoms info energy_per_atoms

      boundary = ExTinyMD.Q2dBoundary(L, L, L)

      ratio = (n_atoms / 1000)^(1 / 3)
      N_real = Int.(ceil.(ratio .* N_real0))
      extra_pad_ratio = Int(ceil(extra_pad_ratio_intial * ratio))

      M_mid = proper_M(eta, L, uspara)
      interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Qs, 0.0, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Ql, Rz_0, Q_0; preset = preset, ϵ = 1.0)
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