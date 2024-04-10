default: init update

init:
	julia --project=. -e 'using Pkg; Pkg.add(url = "https://github.com/HPMolSim/FastSpecSoG.jl"); Pkg.add(url = "https://github.com/HPMolSim/EwaldSummations.jl.git");  Pkg.instantiate(); Pkg.precompile();'

update:
	julia --project=. -e 'using Pkg; Pkg.update(); Pkg.precompile();'

run_T0:
	rm data/Acc_T0_short.csv
	touch data/Acc_T0_short.csv
	echo "preset,M_mid,total,near,error_rel" >> data/Acc_T0_short.csv
	julia --project=. -t 256 data_scripts/T0_Short.jl

run_T1:
	rm data/Acc_T1_Kspace.csv
	touch data/Acc_T1_Kspace.csv
	echo "preset,k_max,energy_exact,energy_k,error,error_rel" >> data/Acc_T1_Kspace.csv
	julia --project=. -t 256 data_scripts/T1_1_Kspace.jl

run_T2_1:
	rm data/Acc_T2_1_Mmid.csv
	touch data/Acc_T2_1_Mmid.csv
	echo "extra_pad_ratio,M_mid,error,error_rel" >> data/Acc_T2_1_Mmid.csv
	julia --project=. data_scripts/T2_1_Mmid.jl

run_T2_2:
	rm data/Acc_T2_2_Extrapadding.csv
	touch data/Acc_T2_2_Extrapadding.csv
	echo "extra_pad_ratio,M_mid,error,error_rel" >> data/Acc_T2_2_Extrapadding.csv
	julia --project=. data_scripts/T2_2_Extrapadding.jl

run_T3_1:
	rm  data/Acc_T3_1_WinFunc_ES.csv
	touch data/Acc_T3_1_WinFunc_ES.csv
	echo "Window_type,energy_exact,error,error_rel,w" >> data/Acc_T3_1_WinFunc_ES.csv
	rm data/Acc_T3_1_WinFunc_Gaussian.csv
	touch data/Acc_T3_1_WinFunc_Gaussian.csv
	echo "Window_type,energy_exact,error,error_rel,w" >> data/Acc_T3_1_WinFunc_Gaussian.csv
	rm data/Acc_T3_1_WinFunc_PKB.csv
	touch data/Acc_T3_1_WinFunc_PKB.csv
	echo "Window_type,energy_exact,error,error_rel,w" >> data/Acc_T3_1_WinFunc_PKB.csv
	julia --project=. data_scripts/T3_1_WinFunc_ES.jl
	julia --project=. data_scripts/T3_1_WinFunc_Gaussian.jl
	julia --project=. data_scripts/T3_1_WinFunc_PKB.jl

run_T3_2:
	rm data/Acc_T3_2_Nspec_PKB.csv
	touch data/Acc_T3_2_Nspec_PKB.csv
	echo "Window_type,energy_exact,error,error_rel,N_real" >> data/Acc_T3_2_Nspec_PKB.csv
	rm data/Acc_T3_2_Nspec_Gaussian.csv
	touch data/Acc_T3_2_Nspec_Gaussian.csv
	echo "Window_type,energy_exact,error,error_rel,N_real" >> data/Acc_T3_2_Nspec_Gaussian.csv
	rm data/Acc_T3_2_Nspec_ES.csv
	touch data/Acc_T3_2_Nspec_ES.csv
	echo "Window_type,energy_exact,error,error_rel,N_real" >> data/Acc_T3_2_Nspec_ES.csv
	julia --project=. data_scripts/T3_2_Nspec_PKB.jl
	julia --project=. data_scripts/T3_2_Nspec_Gaussian.jl
	julia --project=. data_scripts/T3_2_Nspec_ES.jl

run_T4_1:
	rm data/Acc_T4_cube_Nxy.csv
	touch data/Acc_T4_cube_Nxy.csv
	echo "Nxy,M_mid,energy_exact,energy_ffct,error,error_rel" >> data/Acc_T4_cube_Nxy.csv
	julia --project=. -t 32 data_scripts/T4_cube_Nxy.jl

run_T4_2:
	rm data/Acc_T4_cube_Rz.csv
	touch data/Acc_T4_cube_Rz.csv
	echo "Rz,M_mid,energy_exact,energy_ffct,error,error_rel" >> data/Acc_T4_cube_Rz.csv
	julia --project=. -t 32 data_scripts/T4_cube_Rz.jl

run_T4_zeroth:
	rm data/Acc_T4_cube_zeroth.csv
	touch data/Acc_T4_cube_zeroth.csv
	echo "Q_0,M_mid,energy_exact,energy_fgt,error_rel" >> data/Acc_T4_cube_zeroth.csv
	julia --project=. -t 256 data_scripts/T4_cube_zeroth.jl

run_T5_Nxy:
	rm data/Acc_T5_thin_Nxy.csv
	touch data/Acc_T5_thin_Nxy.csv
	echo "Nxy,gamma,E_exact,E_thin,error_rel" >> data/Acc_T5_thin_Nxy.csv
	julia --project=. data_scripts/T5_thin_Nxy.jl

run_T5_Rz:
	rm data/Acc_T5_thin_Rz.csv
	touch data/Acc_T5_thin_Rz.csv
	echo "R_z,gamma,E_exact,E_thin,error_rel" >> data/Acc_T5_thin_Rz.csv
	julia --project=. data_scripts/T5_thin_Rz.jl

run_T5_TaylorQ:
	rm data/Acc_T5_thin_TaylorQ.csv
	touch data/Acc_T5_thin_TaylorQ.csv
	echo "TaylorQ,gamma,E_exact,E_thin,error_rel" >> data/Acc_T5_thin_TaylorQ.csv
	julia --project=. data_scripts/T5_thin_TaylorQ.jl

run_T5_w:	
	rm data/Acc_T5_thin_w.csv
	touch data/Acc_T5_thin_w.csv
	echo "w,gamma,E_exact,E_thin,error_rel" >> data/Acc_T5_thin_w.csv
	julia --project=. data_scripts/T5_thin_w.jl

run_T5_zeroth:
	rm data/Acc_T5_thin_zeroth.csv
	touch data/Acc_T5_thin_zeroth.csv
	echo "Q_0,gamma,E_exact,E_fgt,error_rel" >> data/Acc_T5_thin_zeroth.csv
	julia --project=. -t 256 data_scripts/T5_thin_zeroth.jl

run_Ewald2D_cube:
	julia --project=. -t 256 reference/cube_results.jl

run_Ewald2D_thin:
	julia --project=. -t 256 reference/thin_results.jl

run_cube_total_2:
	rm data/Acc_T6_cube_total_2.csv
	touch data/Acc_T6_cube_total_2.csv
	echo "n_atoms,energy_ewald,energy_fssog,error,error_rel" >> data/Acc_T6_cube_total_2.csv
	julia --project=. data_scripts/T6_cube_total_2.jl

run_cube_total_4:
	rm data/Acc_T6_cube_total_4.csv
	touch data/Acc_T6_cube_total_4.csv
	echo "n_atoms,energy_ewald,energy_fssog,error,error_rel" >> data/Acc_T6_cube_total_4.csv
	julia --project=. data_scripts/T6_cube_total_4.jl

run_cube_total_6:
	rm data/Acc_T6_cube_total_6.csv
	touch data/Acc_T6_cube_total_6.csv
	echo "n_atoms,energy_ewald,energy_fssog,error,error_rel" >> data/Acc_T6_cube_total_6.csv
	julia --project=. data_scripts/T6_cube_total_6.jl

plot:
	julia --project=. figs_scripts/fig_convergence.jl
	julia --project=. figs_scripts/fig_mid_acc.jl
	julia --project=. figs_scripts/fig_long_Nxyz.jl
	julia --project=. figs_scripts/fig_windows.jl
	julia --project=. figs_scripts/fig_thin.jl
	julia --project=. figs_scripts/fig_zeroth.jl
