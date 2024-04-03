#! bin/bash

data_file=data/Acc_T1_Kspace.csv

rm -rf $data_file
touch $data_file
echo "uspara,energy_exact,error,error_rel" >> $data_file

julia --project=. data_scripts/T1_1_Kspace.jl