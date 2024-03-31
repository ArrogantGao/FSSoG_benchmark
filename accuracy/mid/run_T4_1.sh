#! bin/bash

data_file=data/Acc_T4_1_WinFunc_ES.csv
#data_file=data/Acc_T4_1_WinFunc_Gaussian.csv
#data_file=data/Acc_T4_1_WinFunc_PKB.csv

rm -rf $data_file
touch $data_file
echo "Window_type,energy_exact,error,error_rel,w" >> $data_file

julia --project=. mid/T4_1_WinFunc_ES.jl
#julia --project=. mid/T4_1_WinFunc_Gaussian.jl
#julia --project=. mid/T4_1_WinFunc_PKB.jl