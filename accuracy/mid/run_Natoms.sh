#! bin/bash

data_file=data/accuracy_mid_Natoms.csv

touch $data_file
echo "extra_pad_ratio,n_atoms,energy_exact,energy_approx,abs_error,relative_error" >> $data_file

julia --project=. mid/accuracy_mid_Natoms.jl