#! bin/bash

data_file=data_for_manu/Acc_T4_1_wideKspace.csv

rm -rf $data_file
touch $data_file
echo "k_max,energy_exact,error,error_rel,M_mid" >> $data_file

julia --project=. test_for_manu/T4_1_wideKspace.jl