default: init update

init:
	julia --project=. -e 'using Pkg; Pkg.develop(url = "https://github.com/HPMolSim/FastSpecSoG.jl"); Pkg.instantiate(); Pkg.precompile();'

update:
	julia --project=. -e 'using Pkg; Pkg.update(); Pkg.precompile();'

run_Ewald2D_reference:
	julia --project=. -t $(nthreads) reference/cube_results.jl
	julia --project=. -t $(nthreads) reference/thin_results.jl
	julia --project=. -t $(nthreads) reference/thin_high_rho.jl

run_T0:
	julia --project=. -t $(nthreads) data_scripts/T0_Short.jl
	julia --project=. -t $(nthreads) data_scripts/T0_total.jl
	julia --project=. -t $(nthreads) data_scripts/T0_energy.jl
	julia --project=. -t $(nthreads) data_scripts/T0_force.jl

run_T1:
	julia --project=. -t $(nthreads) data_scripts/T1_1_Kspace.jl

run_T2:
	julia --project=. data_scripts/T2_1_Mmid.jl
	julia --project=. data_scripts/T2_2_Extrapadding.jl

run_T3_1:
	julia --project=. data_scripts/T3_1_WinFunc_ES.jl
	julia --project=. data_scripts/T3_1_WinFunc_Gaussian.jl
	julia --project=. data_scripts/T3_1_WinFunc_PKB.jl
	julia --project=. data_scripts/T3_2_Nspec_PKB.jl
	julia --project=. data_scripts/T3_2_Nspec_Gaussian.jl
	julia --project=. data_scripts/T3_2_Nspec_ES.jl

run_T4:
	julia --project=. data_scripts/T4_cube_Nxy.jl
	julia --project=. data_scripts/T4_cube_Rz.jl
	julia --project=. data_scripts/T4_cube_zeroth.jl

run_T5:
	julia --project=. data_scripts/T5_thin_Nxy.jl
	julia --project=. data_scripts/T5_thin_Rz.jl
	julia --project=. data_scripts/T5_thin_TaylorQ.jl
	julia --project=. data_scripts/T5_thin_w.jl
	julia --project=. data_scripts/T5_thin_zeroth.jl

run_cube_total:
	julia --project=. data_scripts/T6_cube_total_2.jl
	julia --project=. data_scripts/T6_cube_total_4.jl
	julia --project=. data_scripts/T6_cube_total_6.jl

run_thin_total:
	julia --project=. data_scripts/T6_thin_total_2.jl
	julia --project=. data_scripts/T6_thin_total_4.jl
	julia --project=. data_scripts/T6_thin_total_6.jl

run_high_dense:
	julia --project=. data_scripts/T7_cube_total_2_highdense.jl
	julia --project=. data_scripts/T7_thin_total_6_highdense.jl

plot:
	julia --project=. figs_scripts/fig_convergence.jl
	julia --project=. figs_scripts/fig_mid_acc.jl
	julia --project=. figs_scripts/fig_long_Nxyz.jl
	julia --project=. figs_scripts/fig_windows.jl
	julia --project=. figs_scripts/fig_thin.jl
	julia --project=. figs_scripts/fig_zeroth.jl
	julia --project=. figs_scripts/fig_total.jl
