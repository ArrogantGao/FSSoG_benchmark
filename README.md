# Repo for benchmarking

## How to init
Type the following commands in bash
```
git clone git@github.com:ArrogantGao/FSSoG_benchmark.git
cd FSSoG_benchmark
julia --project=.
```
Then in Julia REPL
```
(FSSoG_benchmark) pkg> add https://github.com/HPMolSim/FastSpecSoG.jl

(FSSoG_benchmark) pkg> instantiate

(FSSoG_benchmark) pkg> resolve
  No Changes to `~/temp/FSSoG_benchmark/Project.toml`
  No Changes to `~/temp/FSSoG_benchmark/Manifest.toml`

julia> exit()
```
To run the scripts, for example
```
cd FFCT
julia --project=.. FFCT_profile.jl
```
