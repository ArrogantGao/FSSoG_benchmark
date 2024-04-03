#! bin/bash

data_file1=data_for_manu/Acc_T3_1_WinFunc_ES.csv
data_file2=data_for_manu/Acc_T3_1_WinFunc_Gaussian.csv
data_file3=data_for_manu/Acc_T3_1_WinFunc_PKB.csv

rm -rf $data_file1
touch $data_file1
echo "Window_type,energy_exact,error,error_rel,w" >> $data_file1

rm -rf $data_file2
touch $data_file2
echo "Window_type,energy_exact,error,error_rel,w" >> $data_file2

rm -rf $data_file3
touch $data_file3
echo "Window_type,energy_exact,error,error_rel,w" >> $data_file3


julia --project=. test_for_manu/T3_1_WinFunc_ES.jl
julia --project=. test_for_manu/T3_1_WinFunc_Gaussian.jl
julia --project=. test_for_manu/T3_1_WinFunc_PKB.jl