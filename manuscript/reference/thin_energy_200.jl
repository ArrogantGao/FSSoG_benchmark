using EwaldSummations, ExTinyMD, JLD2, FastSpecSoG
using Random
Random.seed!(123)

function main()

    n_atoms0 = 1000
    Lx0 = 30.0
    Lz = 0.3

    s = 6.0
    ratios = [10^i for i in 2.5:0.5:3.0]

    for ratio in ratios
        n_atoms = Int(ceil(n_atoms0 * ratio))
        if isodd(n_atoms)
            n_atoms += 1
        end
    
        Lx = Lx0 * ratio^(1/2)
        
        boundary = ExTinyMD.Q2dBoundary(Lx, Lx, Lz)
    
        atoms = Vector{Atom{Float64}}()
        for i in 1:n_atoms÷2
            push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
        end
    
        for i in n_atoms÷2 + 1 : n_atoms
            push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
        end
    
        info = SimulationInfo(n_atoms, atoms, (0.0, Lx, 0.0, Lx, 0.0, Lz), boundary; min_r = 0.0, temp = 1.0)
    
        α = s / (Lx / 2 - 0.1)
        Ewald2D_interaction = Ewald2DInteraction(n_atoms, s, α, (Lx, Lx, Lz), ϵ = 1.0)
    
        r_c = Ewald2D_interaction.r_c
    
        @show ratio, n_atoms, Lx, Lz, n_atoms/(Lx^2 * Lz), α, r_c

        p = [p_info.position for p_info in info.particle_info]
        q = [atoms[p_info.id].charge for p_info in info.particle_info]

        @time Es = Ewald2D_short_energy_N(200, Ewald2D_interaction, p, q) 
        @time El = Ewald2D_long_energy_N(200, Ewald2D_interaction, p, q)

        energy_per_atoms = Es + El
    
        @save "reference/thin_200/n_$(n_atoms).jld2" n_atoms Lx Lz atoms info energy_per_atoms
    end

end

main()