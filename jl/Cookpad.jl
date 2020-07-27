using Pkg
Pkg.activate(".")
using SimpleHypergraphs, Random, ProgressBars, Plots, JLD2
include("HypergraphClustering.jl")
include("ClusteringUtil.jl")


function build_trimcookpad(fname)
    io = open(fname, "r")
    lnum = 0
    recipe = []
    h = Hypergraph{Int}(75874, 54812)
    recipe_dict = Dict{Int,Int}()
    for line in tqdm(eachline(io)) #for each he
        lnum += 1
        if(lnum == 1)
            continue
        end

        line = split(line, " ")

        ingredient_list = split.(line[2:end-2], ":")

        for i in ingredient_list
            ingredient = parse(Int, String(i[1]))
            h[ingredient, lnum-1] = 1
        end
    end
    close(io)

    return h

end

@time const trim_cookpad = build_trimcookpad("../ingredients_trim15.blbl")

@time ms, ph, eph, bcn, uf = he_clustering(trim_cookpad, 1, my_mod, freq=100)
@save "cookpad_heclus.jld2" ms ph eph bcn uf
