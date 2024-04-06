#! bin/bash

data_file=data/Acc_T1_Kspace.csv

rm -rf $data_file
touch $data_file
echo "preset,k_max,energy_exact,energy_k,error,error_rel" >> $data_file

julia --project=. -p 127 -t 256 data_scripts/T1_1_Kspace.jl