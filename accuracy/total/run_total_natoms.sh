#! bin/bash

data_file=data/accuracy_total_natoms.csv

rm $data_file
touch $data_file
echo "n_atoms,E_exact,E_fssog,abs_error,relative_error" >> $data_file

julia --project=. total/accuracy_total_natoms.jl