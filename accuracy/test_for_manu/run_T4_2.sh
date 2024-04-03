#! bin/bash

data_file=data_for_manu/Acc_T4_2_Cheb.csv

rm -rf $data_file
touch $data_file
echo "N_z,M_mid,E_exact,E_FFCT,abs_error,relative_error" >> $data_file

julia --project=. test_for_manu/T4_2_Cheb.jl