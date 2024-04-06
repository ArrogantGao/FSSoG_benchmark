#! bin/bash

data_file=data/Acc_T4_1_wideKspace.csv

rm -rf $data_file
touch $data_file
echo "k_max,energy_exact,error,error_rel,M_mid" >> $data_file

julia --project=. -p 127 -t 256 data_scripts/T4_1_wideKspace.jl