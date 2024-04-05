#! bin/bash

data_file=data/Acc_T5_thin_w.csv

rm -rf $data_file
touch $data_file
echo "w,error_rel_pkb,error_ES,error_Gaussian" >> $data_file

julia --project=. data_scripts/T5_thin_w.jl