#!/bin/bash

output_file_cube=$1
output_file_thin=$2
output_file_cube_highdense=$3
output_file_thin_highdense=$4

touch "$output_file_cube"
echo 'preset,n_atoms,t_short,t_mid,t_long,t_total' >> $output_file_cube

touch "$output_file_cube_highdense"
echo 'preset,n_atoms,t_short,t_mid,t_long,t_total' >> $output_file_cube_highdense

touch "$output_file_thin"
echo 'preset,n_atoms,t_short,t_long,t_total' >> $output_file_thin

touch "$output_file_thin_highdense"
echo 'preset,n_atoms,t_short,t_long,t_total' >> $output_file_thin_highdense

# Start recording the benchmark execution time
start_time=$(date +%s)

echo -e "Starting benchmark..."
# julia --threads 1 --project=@. run_cube_6_highdense.jl --output-file "$output_file_cube_highdense"
julia --threads 1 --project=@. run_thin_6_highdense.jl --output-file "$output_file_thin_highdense"
# julia --threads 1 --project=@. run_cube_2.jl --output-file "$output_file_cube"
# julia --threads 1 --project=@. run_cube_4.jl --output-file "$output_file_cube"
# julia --threads 1 --project=@. run_cube_6.jl --output-file "$output_file_cube"
# julia --threads 1 --project=@. run_thin_2.jl --output-file "$output_file_thin"
# julia --threads 1 --project=@. run_thin_4.jl --output-file "$output_file_thin"
# julia --threads 1 --project=@. run_thin_6.jl --output-file "$output_file_thin"
echo -e "Benchmark finished"

# Log the benchmark execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))
formatted_time=$(date -u -d @"$execution_time" +%H:%M:%S)
echo "Benchmark execution time: $formatted_time"

echo -e "Generating plots..."
julia --project=@. generate-graph-cube.jl --data-file $output_file_cube
julia --project=@. generate-graph-thin.jl --data-file $output_file_thin
echo -e "Plot generated"
