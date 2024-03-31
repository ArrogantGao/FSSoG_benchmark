#! bin/bash

data_file=data/Acc_T1_2_wideKspace.csv

rm -rf $data_file
touch $data_file
echo "k_max,energy_exact,error,error_rel,M_mid" >> $data_file

julia --project=. long/T1_2_wideKspace.jl