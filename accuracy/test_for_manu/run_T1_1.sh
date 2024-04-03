#! bin/bash

data_file=data_for_manu/Acc_T1_Kspace.csv

rm -rf $data_file
touch $data_file
echo "uspara,energy_exact,error,error_rel" >> $data_file

julia --project=. test_for_manu/T1_1_Kspace.jl