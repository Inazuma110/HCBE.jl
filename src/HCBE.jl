module HCBE
using SimpleHypergraphs, Random,  ProgressBars, LightGraphs, ProgressMeter, GraphPlot, StatsBase, JSON, StatsPlots

include("./algorithms/HypergraphClustering.jl")
include("./algorithms/Weighted.jl")
include("./algorithms/ClusteringUtil.jl")
include("./algorithms/CreateHypergraph.jl")
include("./IO.jl")
end # module
