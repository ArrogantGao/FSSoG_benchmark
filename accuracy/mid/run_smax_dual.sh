#! bin/bash

data_file=data/accuracy_mid_smax_dual.csv

touch $data_file
echo "extra_pad_ratio,M_mid,error,error_rel" >> $data_file

julia --project=. mid/accuracy_mid_smax_dual.jl