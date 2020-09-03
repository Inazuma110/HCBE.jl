using Pkg
using Random, SimpleHypergraphs, StatsBase, ProgressMeter, ProgressBars

# clusters cluster num
# npcs node per cluster 1つあたりの頂点数
function create_hypergraph(npcs, hepcs, he_rate=0.5, noise_rate=0.1)
  node_num = 0
  he_num = 0
  clusters = length(npcs)
  training_data = []
  for npc in npcs
    cluster = Set([node_num+i for i in 1:npc])
    training_data = vcat(training_data, cluster)
    node_num += npc
  end

  for hepc in hepcs
    he_num += hepc
  end
  h = Hypergraph{Int64}(node_num, he_num)
  step = 1
  fin_node = 1
  is_noises = [[false for j in 1:node_num+he_num] for i in 1:node_num+he_num]

  for (npc, hepc) in zip(npcs, hepcs)
    for i in 1:hepc
      # he_size = rand(1:npc)
      he = randsubseq(fin_node:fin_node+npc-1, he_rate)
      he = [(i in he ? 1 : nothing) for i in 1:nhv(h)]
      for (j, v) in enumerate(he)
        h[j, step] = v
      end
      step += 1
    end
    fin_node += npc
  end

  for node in 1:nhv(h)
    is_noise = rand(1)[1] <= noise_rate
    # println(rand(1)[1])
    if !is_noise continue end

    included_cluster_num = StatsBase.sample(1:length(npcs))
    included_clusters = randperm(included_cluster_num)

    for cluster in included_clusters
      included_he = randsubseq(1:nhe(h), rand(1)[1])
      for he in included_he
        if hepcs[cluster] * (cluster-1) >= he continue end
        if hepcs[cluster] * cluster < he continue end
        h[node, he] = 1
      end
      # h[node, cluster*npcs[1] .+ included_he] .= 1
      # is_noises[node][node_num.+included_he] .= true
    end

  end

  return h, is_noises, training_data
end

function build_trimcookpad(fname="../ingredients_trim15.blbl")
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


function build_youtube(fname="../youtube_giant.blbl")
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

function build_wiki(fname="../wiki/enwiki-2013.txt")
    io = open(fname, "r")
    lnum = 0
    recipe = []
    h = Hypergraph{Int}(4203323, 4203323)
    recipe_dict = Dict{Int,Int}()
    for line in tqdm(eachline(io))
      lnum += 1
      if(lnum <= 4)
        continue
      end

      line = parse.(Int, split(line, " "))
      h[line[1]+1, line[2]+1] = 1
    end
    close(io)
    return h
end

function build_amazon(fname="../amazon-meta.txt")
    io = open(fname, "r")
    lnum = 0
    id_s = Set()
    customer_s = Set()
    d = Dict()
    data = Dict()
    delimitor = r",|\t|:| "
    for line in tqdm(eachline(io))
      lnum += 1
      if lnum <= 7 continue end

      line = split(line, delimitor, keepempty=false)
      if isempty(line)
        if haskey(d, "Id") push!(id_s, d["Id"]) end
        for (k, v) in d
          if '-' in k && length(v) >= 2 push!(customer_s, v[2]) end
        end
      else d[line[1]] = line[2:end] end
    end
    close(io)
    return id_s, customer_s
end

function build_dblp(fname="../com-dblp.all.cmty.txt")
  io = open(fname, "r")
  lnum = 0
  d = Dict([])
  s = Set()
  for line in tqdm(eachline(io))
    line = parse.(Int, split(line, "\t"))
    for node in line
      push!(s, node)
    end
  end
  close(io)

  d = Dict([node=>i for (i, node) in enumerate(s)])

  io = open(fname, "r")
  h = Hypergraph{Int64}(length(s), 13477)
  for line in tqdm(eachline(io))
    lnum += 1

    line = parse.(Int, split(line, "\t"))
    for node in line
      h[d[node], lnum] = 1
    end
  end
  close(io)
  return h
end

function build_sample()
  h1 = Hypergraph{Int64}(6,9)
  h1[1:3,1] .= 1
  h1[1:2,2] .= 1
  h1[2:3,3] .= 1
  h1[1, 4] = 1
  h1[3, 4] = 1
  h1[3:4, 5] .= 1
  h1[4:6,6] .= 1
  h1[4:5,7] .= 1
  h1[5:6, 8] .= 1
  h1[4, 9] = 1
  h1[6, 9] = 1
  # h1[3, 5:9] .= 1
  return h1
end

