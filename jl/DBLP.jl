using Pkg
Pkg.activate(".")
using SimpleHypergraphs, Random, ProgressBars, Plots, JLD2
include("HypergraphClustering.jl")
include("ClusteringUtil.jl")
include("./CreateHypergraph.jl")


@time const dblp = build_dblp()

@time ms, ph, bp, ufh = clustering3(dblp, 1, my_mod, freq=10000)
@save "cookpad_heclus.jld2" ms ph bp ufh
