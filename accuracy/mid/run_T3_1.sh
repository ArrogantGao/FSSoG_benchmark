#! bin/bash

data_file=data/Acc_T3_1_window.csv

rm -rf $data_file
touch $data_file
echo "extra_pad_ratio,N_atoms,E_exact,E_approx,error,error_rel" >> $data_file

julia --project=. mid/T3_1_window.jl