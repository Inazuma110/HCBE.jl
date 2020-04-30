using Pkg
Pkg.activate(".")
using SimpleHypergraphs, Random, ProgressBars, Plots, JLD2
include("HypergraphClustering.jl")

function build_youtube(fname)
    io = open(fname, "r")
    lnum = 0
    recipe = []
    h = Hypergraph{Int}(45352, 13251)
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

@time const youtube = build_youtube("../youtube_giant.blbl")
# @time uf, ms, cs = clustering4(youtube)
cfm = CFModularityCNMLike(10000)
@time bm, bp, mhist = findcommunities(youtube, cfm)

@save "youtube_jl.dat" bm bp mhist
