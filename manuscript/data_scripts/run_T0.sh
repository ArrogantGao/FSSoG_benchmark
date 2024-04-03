#! bin/bash

data_file=data/Acc_T0_short.csv

rm -rf $data_file
touch $data_file
echo "uspara, error" >> $data_file

julia --project=. data_scripts/T0_Short.jl