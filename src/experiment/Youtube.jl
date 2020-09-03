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

function find_first(c::Array{Set{Int}}, vals)
    for i in 1:length(c)
        for v in vals
            v in c[i] && return i
        end
    end
    throw("None of values in $vals found")
end

function f(h::Hypergraph, method::CFModularityCNMLike; mod_f=my_mod)
    # ha = HypergraphAggs(h)
    best_modularity = 0
    comms = [Set(i) for i in 1:nhv(h)]
    mod_history = Vector{Float64}(undef, method.reps)
    @showprogress for rep in 1:method.reps
        he = rand(1:nhe(h))
        vers = collect(keys(getvertices(h, he)))
        if length(vers) == 0
            continue
        end;
        c = deepcopy(comms)
        i0 = find_first(c, vers)
        max_i = length(c)
        i_cur = i0
        while i_cur < max_i
            i_cur += 1
            if length(intersect(c[i_cur],vers)) > 0
                union!(c[i0], c[i_cur])
                c[i_cur]=c[max_i]
                max_i += -1
            end
        end
        resize!(c,max_i)
        m = mod_f(h, c)
        if m > best_modularity
            best_modularity = m
            comms = c
        end
        mod_history[rep] = best_modularity
    end
    return (bm=best_modularity, bp=comms, mod_history=mod_history)
end

@time const youtube = build_youtube("../youtube_giant.blbl")
# @time uf, ms, cs = clustering4(youtube)
cfm = CFModularityCNMLike(500)
@time bm, bp, mhist = f(youtube, cfm)

@save "youtube_jl.dat" bm bp mhist
