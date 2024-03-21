using FastSpecSoG, ExTinyMD, BenchmarkTools, Plots, CSV, DataFrames, ArgParse

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

    n_atoms0 = 100
    L0 = 100.0
    N_real0 = (64, 64, 64)


    w = (16, 16, 16)
    β = 5.0 .* w
    extra_pad_ratio = 2
    cheb_order = 10
    preset = 2
    M_mid = 6
    N_grid = (4, 4, 16)
    Q = 10
    r_c = 10.0

    for ratio in [2^i for i in 0:0.5:2]

        n_atoms = Int(ceil(n_atoms0 * ratio^3))
        if isodd(n_atoms)
            n_atoms += 1
        end
        L = L0 * ratio

        boundary = ExTinyMD.Q2dBoundary(L, L, L)

        atoms = Vector{Atom{Float64}}()
        for i in 1:n_atoms÷2
            push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
        end

        for i in n_atoms÷2 + 1 : n_atoms
            push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
        end

        info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

        N_real = Int.(ceil.(ratio .* N_real0))

        @info "FSSoG, n_atoms = $n_atoms"
        interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Q, 0.5, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Q, SoePara16(); preset = preset, ϵ = 1.0)
        neighbor = CellList3D(info, r_c, boundary, 1)
        
        for i in 1:interaction.n_atoms
            interaction.charge[i] = atoms[info.particle_info[i].id].charge
            interaction.position[i] = info.particle_info[i].position.coo
        end

        t_short = @belapsed energy_short($interaction, $neighbor)
        t_mid = @belapsed energy_mid($interaction)
        t_long = @belapsed energy_long($interaction)

        @show n_atoms, t_short, t_mid, t_long

        df = DataFrame(n_atoms = n_atoms, t_short = t_short, t_mid = t_mid, t_long = t_long)
        CSV.write(output_file, df, append = true)
    end
end

main()