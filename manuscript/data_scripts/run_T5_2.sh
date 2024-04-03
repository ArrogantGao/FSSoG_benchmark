#! bin/bash

data_file=data/Acc_T5_2_ThinCheb.csv

rm -rf $data_file
touch $data_file
echo "preset,Rz,E_thin_k,E_direct_k,error_rel" >> $data_file

julia --project=. data_scripts/T5_2_ThinCheb.jl