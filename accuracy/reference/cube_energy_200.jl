using EwaldSummations, ExTinyMD, JLD2
using Random
Random.seed!(123)

function main()

    n_atoms0 = 1000
    L0 = 20.0

    s = 6.0
    ratios = [10.0^i for i in 0.0:0.5:3]

    for ratio in ratios
        n_atoms = Int(ceil(n_atoms0 * ratio))
        if isodd(n_atoms)
            n_atoms += 1
        end

        L = L0 * ratio^(1/3)
        
        boundary = ExTinyMD.Q2dBoundary(L, L, L)

        atoms = Vector{Atom{Float64}}()
        for i in 1:n_atoms÷2
            push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
        end

        for i in n_atoms÷2 + 1 : n_atoms
            push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
        end

        @info "system init"
        info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 0.0, temp = 1.0)
        @info "init done"

        #s / α = L / 2 - 0.1
        α = s / (L / 2 - 0.1)
        Ewald2D_interaction = Ewald2DInteraction(n_atoms, s, α, (L, L, L), ϵ = 1.0)

        r_c = Ewald2D_interaction.r_c

        @show ratio, n_atoms, L, n_atoms/L^3, α, r_c

        # @info "finding neighbors"
        # Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
        # @info "finding neighbors finished"

        p = [p_info.position for p_info in info.particle_info]
        q = [atoms[p_info.id].charge for p_info in info.particle_info]

        @time Es = Ewald2D_short_energy_N(200, Ewald2D_interaction, p, q) 
        @time El = Ewald2D_long_energy_N(200, Ewald2D_interaction, p, q)

        energy_per_atoms = Es + El

        @save "reference/cube_200/n_$(n_atoms).jld2" n_atoms L atoms info energy_per_atoms
    end
end

main()