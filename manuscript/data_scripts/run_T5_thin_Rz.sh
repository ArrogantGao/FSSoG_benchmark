#! bin/bash

data_file=data/Acc_T5_thin_Rz.csv

rm -rf $data_file
touch $data_file
echo "R_z,E_exact,E_thin,error_rel" >> $data_file

julia --project=. data_scripts/T5_thin_Rz.jl