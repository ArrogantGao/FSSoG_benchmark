#! bin/bash

data_file=data/Acc_T1_2_WinFunc_ES.csv

rm -rf $data_file
touch $data_file
echo "Window_type,energy_exact,error,error_rel,w" >> $data_file

julia --project=. mid/T1_2_WinFunc_ES.jl