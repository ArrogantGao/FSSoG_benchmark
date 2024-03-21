#! bin/bash

data_file=data/accuracy_long_Nz.csv

touch $data_file
echo "N_z,M_mid,E_exact,E_FFCT,abs_error,relative_error" >> $data_file

julia --project=. long/accuracy_long_Nz.jl