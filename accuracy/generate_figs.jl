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