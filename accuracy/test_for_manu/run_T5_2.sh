#! bin/bash

data_file=data_for_manu/Acc_T5_2_ThinCheb.csv

rm -rf $data_file
touch $data_file
echo "preset,Rz,E_thin_k,E_direct_k,error_rel" >> $data_file

julia --project=. test_for_manu/T5_2_ThinCheb.jl