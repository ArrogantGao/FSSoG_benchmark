#! bin/bash

#data_file=data/Acc_T4_2_Nspec_PKB.csv
data_file=data/Acc_T4_2_Nspec_Gaussian.csv
#data_file=data/Acc_T4_2_Nspec_ES.csv

rm -rf $data_file
touch $data_file
echo "Window_type,energy_exact,error,error_rel,N_real" >> $data_file

#julia --project=. mid/T4_2_Nspec_PKB.jl
julia --project=. mid/T4_2_Nspec_Gaussian.jl
#julia --project=. mid/T4_2_Nspec_ES.jl