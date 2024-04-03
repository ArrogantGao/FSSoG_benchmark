#! bin/bash

data_file1=data_for_manu/Acc_T3_2_Nspec_PKB.csv
data_file2=data_for_manu/Acc_T3_2_Nspec_Gaussian.csv
data_file3=data_for_manu/Acc_T3_2_Nspec_ES.csv

rm -rf $data_file1
touch $data_file1
echo "Window_type,energy_exact,error,error_rel,N_real" >> $data_file1

rm -rf $data_file2
touch $data_file2
echo "Window_type,energy_exact,error,error_rel,N_real" >> $data_file2

rm -rf $data_file3
touch $data_file3
echo "Window_type,energy_exact,error,error_rel,N_real" >> $data_file3

julia --project=. test_for_manu/T3_2_Nspec_PKB.jl
julia --project=. test_for_manu/T3_2_Nspec_Gaussian.jl
julia --project=. test_for_manu/T3_2_Nspec_ES.jl