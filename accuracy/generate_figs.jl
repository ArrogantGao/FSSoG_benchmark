using CSV, DataFrames, Plots, LaTeXStrings

function mid_error_M_mid()
    df = CSV.read("data/accuracy_mid_smax_w16.csv", DataFrame)

    grouped_dfs = groupby(df, :extra_pad_ratio)

    # change the df into five arrays
    extra_pad_ratio = unique(df.extra_pad_ratio)
    M_mid = unique(df.M_mid)
    errors = [dfi[:, 3] for dfi in grouped_dfs]

    plot(xlabel = L"M_{mid}", ylabel = "Error", title = L"Accuracy of mid-range part to $M_{mid}$", dpi = 500, xlim = [0, 20.5], ylim = [-15, 1.0])
    for i in 1:length(extra_pad_ratio)
        plot!(M_mid, log10.(errors[i]), label = "λ = $(extra_pad_ratio[i])", marker = :circle)
    end
    plot!(legend = :bottomright)
    savefig("figs/accuracy_mid_Mmid.png")
end

mid_error_M_mid()

function mid_Natoms()
    df = CSV.read("data/accuracy_mid_Natoms.csv", DataFrame)
    grouped_dfs = groupby(df, :extra_pad_ratio)
    extra_pad_ratio = unique(df.extra_pad_ratio)
    N_atoms = unique(df.n_atoms)

    relative_errors = [dfi[:, 6] for dfi in grouped_dfs]

    plot(xlabel = "log10(N_atoms)", ylabel = "log10(Relative Error)", title = "Accuracy of mid-range part to N_atoms", dpi = 500, xlim = [1.8, 4], ylim = [-15, 1.0])
    for i in 1:length(extra_pad_ratio)
        plot!(log10.(N_atoms), log10.(relative_errors[i]), label = "λ = $(extra_pad_ratio[i])", marker = :circle)
    end
    plot!(legend = :topright)
    savefig("figs/accuracy_mid_Natoms.png")
end

mid_Natoms()

function long_Nz()
    df0 = CSV.read("data/accuracy_long_Nz.csv", DataFrame)
    M_mid_array = unique(df0.M_mid)
    plot(xlabel = "N_z", ylabel = "log10(Relative error)", legend = :topright, dpi = 500)
    for i in 1:2:length(M_mid_array)
        M_mid = M_mid_array[i]
        df = filter(row -> row.M_mid == M_mid, df0)
        plot!(df.N_z, log10.(df.relative_error), label = "M_mid = $M_mid", marker = :circle)
    end
    savefig("figs/accuracy_long_Nz.png")
end

long_Nz()

function long_Mmid()
    df = CSV.read("data/accuracy_long_Mmid.csv", DataFrame)
    plot(xlabel = L"M_{mid}", ylabel = "contribution", title = L"Contribution by $M_{mid}$", dpi = 500, xlim = [0, 30.5], ylim = [-15, 1.0])
    for i in 1:5
        dfi = filter(row -> row.preset == i, df)
        plot!(dfi.M_mid, log10.(abs.(dfi.E_exact)), label = "preset = $(i)", marker = :circle)
    end
    plot!(legend = :bottomright)
    savefig("figs/accuracy_long_Mmid_k.png")

    plot(xlabel = L"M_{mid}", ylabel = "contribution", title = L"Contribution by $M_{mid}$, zeroth order", dpi = 500, xlim = [0, 150], ylim = [-15, 1.0])
    for i in 1:5
        dfi = filter(row -> row.preset == i, df)
        plot!(dfi.M_mid, log10.(abs.(dfi.E_0)), label = "preset = $(i)", marker = :circle)
    end
    plot!(legend = :topright)
    savefig("figs/accuracy_long_Mmid_0.png")
end

long_Mmid()

function long_Nxy()
    df0 = CSV.read("data/accuracy_long_Nxy.csv", DataFrame)
    M_mid_array = unique(df0.M_mid)
    plot(xlabel = "log2(N_xy)", ylabel = "log10(Relative error)", legend = :topright, dpi = 500)
    for i in 1:length(M_mid_array)
        M_mid = M_mid_array[i]
        df = filter(row -> row.M_mid == M_mid, df0)
        plot!(log2.(df.N_xy), log10.(df.relative_error), label = "M_mid = $M_mid", marker = :circle)
    end
    savefig("figs/accuracy_long_Nxy.png")
end

long_Nxy()

function total_natoms()
    df = CSV.read("data/accuracy_total_natoms_sog.csv", DataFrame)
    plot(xlabel = "log10(n_atoms)", ylabel = "log10(Relative error)", legend = :topright, dpi = 500, ylim = [-15, 1.0])
    plot!(log10.(df.n_atoms), log10.(df.relative_error), marker = :circle, label = "preset = 2")
    savefig("figs/accuracy_total_natoms.png")
end

total_natoms()