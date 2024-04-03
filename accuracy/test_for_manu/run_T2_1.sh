#! bin/bash

data_file=data_for_manu/Acc_T2_1_Mmid.csv

rm -rf $data_file
touch $data_file
echo "extra_pad_ratio,M_mid,error,error_rel" >> $data_file

julia --project=. test_for_manu/T2_1_Mmid.jl