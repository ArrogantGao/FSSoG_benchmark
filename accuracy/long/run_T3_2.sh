#! bin/bash

data_file=data/Acc_T3_2_Cheb_additional.csv

rm -rf $data_file
touch $data_file
echo "N_z,M_mid,E_exact,E_FFCT,abs_error,relative_error" >> $data_file

julia --project=. long/T3_2_Cheb.jl