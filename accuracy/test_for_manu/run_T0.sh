#! bin/bash

data_file=data_for_manu/Acc_T0_short.csv

rm -rf $data_file
touch $data_file
echo "uspara, error" >> $data_file

julia --project=. test_for_manu/T0_Short.jl