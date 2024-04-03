#! bin/bash

data_file=data/Acc_T2_2_Extrapadding.csv

rm -rf $data_file
touch $data_file
echo "extra_pad_ratio,M_mid,error,error_rel" >> $data_file

julia --project=. data_scripts/T2_2_Extrapadding.jl