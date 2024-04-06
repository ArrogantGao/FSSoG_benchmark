#! bin/bash

data_file=data/Acc_T4_cube_Rz.csv

rm -rf $data_file
touch $data_file
echo "Rz,M_mid,energy_exact,energy_ffct,error,error_rel" >> $data_file

julia --project=. -p 32 data_scripts/T4_cube_Rz.jl