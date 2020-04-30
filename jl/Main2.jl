using Pkg
Pkg.activate(".")
using SimpleHypergraphs, Random, ProgressBars, Plots, JLD2
include("HypergraphClustering.jl")

function build_cookpad(fname)
    io = open(fname, "r")
    lnum = 0
    recipe = []
    h = Hypergraph{Int}(487568, 1712897)
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

@time const cookpad = build_cookpad("../ingredients_giant.blbl")

cfm = CFModularityCNMLike(1)
@time bm, bp, mod_history = findcommunities(cookpad, cfm)

@save "res2.dat" bm bp mod_history
