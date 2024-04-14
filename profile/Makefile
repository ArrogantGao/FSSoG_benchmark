default: init update

init:
	julia --project=. -e 'using Pkg; Pkg.add(url = "https://github.com/HPMolSim/FastSpecSoG.jl"); Pkg.add(url = "https://github.com/HPMolSim/EwaldSummations.jl.git");  Pkg.instantiate(); Pkg.precompile();'

update:
	julia --project=. -e 'using Pkg; Pkg.update(); Pkg.precompile();'

profile_thin:
	julia --project=. thin/run_thin_6.jl
	mv statprof thin/prof_results_6