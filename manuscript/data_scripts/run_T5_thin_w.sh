#! bin/bash

data_file=data/Acc_T5_thin_w.csv

rm -rf $data_file
touch $data_file
echo "w,gamma,E_exact,E_thin,error_rel" >> $data_file

julia --project=. data_scripts/T5_thin_w.jl