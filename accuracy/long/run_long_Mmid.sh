#! bin/bash

data_file=data/accuracy_long_Mmid.csv

touch $data_file
echo "preset,M_mid,E_exact,E_0" >> $data_file

julia --project=. long/accuracy_long_Mmid.jl