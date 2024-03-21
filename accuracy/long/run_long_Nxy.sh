#! bin/bash

data_file=data/accuracy_long_Nxy.csv

touch $data_file
echo "N_xy,M_mid,E_exact,E_FFCT,abs_error,relative_error" >> $data_file

julia --project=. long/accuracy_long_Nxy.jl