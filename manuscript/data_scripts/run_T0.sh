#! bin/bash

data_file=data/Acc_T0_short.csv

rm -rf $data_file
touch $data_file
echo "preset,M_mid,total,near,error_rel" >> $data_file

julia --project=. -p 127 -t 256 data_scripts/T0_Short.jl