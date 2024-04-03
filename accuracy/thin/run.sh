#! bin/bash

# data_file=data/thin_long_Mmid.csv
# rm $data_file
# touch $data_file
# echo "M_mid,Nxy,E_thin_k,E_direct_k,error" >> $data_file
# julia --project=. thin/thin_long_Mmid.jl

data_file=data/thin_long_Nxy.csv
rm $data_file
touch $data_file
echo "preset,Nxy,E_thin_k,E_direct_k,error" >> $data_file
julia --project=. thin/thin_long_Nxy.jl

data_file=data/thin_long_Rz.csv
rm $data_file
touch $data_file
echo "preset,Rz,E_thin_k,E_direct_k,error" >> $data_file
julia --project=. thin/thin_long_Rz.jl

data_file=data/thin_long_taylorQ.csv
rm $data_file
touch $data_file
echo "preset,Taylor_Q,E_thin_k,E_direct_k,error" >> $data_file
julia --project=. thin/thin_long_taylorQ.jl
<<<<<<< HEAD
=======

# data_file=data/thin_long_Nxy_loop.csv
# rm $data_file
# touch $data_file
# echo "preset,Nxy,E_thin_k,E_direct_k,error" >> $data_file
# julia --project=. thin/thin_long_Nxy_loop.jl

# data_file=data/thin_long_Nxy_cheb.csv
# rm $data_file
# touch $data_file
# echo "preset,Nxy,E_thin_k,E_direct_k,error" >> $data_file
# julia --project=. thin/thin_long_Nxy_cheb.jl
>>>>>>> main
