#!/bin/bash

output_file=$1

touch "$output_file"
echo 'problem,n_atoms,t_short,t_mid,t_long' >> $output_file

# Start recording the benchmark execution time
start_time=$(date +%s)

echo -e "Starting benchmark..."
julia --threads 1 --project=@. run_energy.jl --output-file "$output_file"
echo -e "Benchmark finished"

# Log the benchmark execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))
formatted_time=$(date -u -d @"$execution_time" +%H:%M:%S)
echo "Benchmark execution time: $formatted_time"

# echo -e "Generating plots..."
# julia --project=@. generate-graph.jl $output_file $complexity_file
# echo -e "Plot generated"