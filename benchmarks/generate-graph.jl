using Plots, CSV, DataFrames, ArgParse

function parse_commandline()

    s = ArgParseSettings()
  
    @add_arg_table s begin
      "--data-file"
      help = "relative file path of the output file"
    end
  
    return parse_args(s)
  
  end

function runtime(data_file)
    df = CSV.read(data_file, DataFrame)
    file_path = joinpath(dirname(data_file), "time_fssog.png")
    n_atoms = df.n_atoms
    t_short = df.t_short
    t_mid = df.t_mid
    t_long = df.t_long
    plot(xlabel = "n_atoms", ylabel = "time (s)", title = "accuracy ~ 1e-4", legend = :topleft, xscale = :log10, yscale = :log10, dpi = 500, xlim = [10^(2.8), 10^(5.2)])
    plot!(n_atoms, 10^(-6) .* n_atoms, label = "linear", linestyle = :dash, color = :black)
    plot!(n_atoms, t_short, label = "short", marker = :circle)
    plot!(n_atoms, t_mid, label = "mid", marker = :circle)
    plot!(n_atoms, t_long, label = "long", marker = :circle)
    @show file_path
    savefig(file_path)
end

data_file = parse_commandline()["data-file"]
runtime(data_file)