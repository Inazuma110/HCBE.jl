using Pkg
using Random, SimpleHypergraphs

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

    # included_num = (rand(Int64, 1)[1] % Int(floor(nhe(h) * noise_rate))) + 1
    # included_he = randperm(nhe(h))[1:included_num]
    included_he = randsubseq(1:nhe(h), he_rate)
    h[node, included_he] .= 1
    is_noises[node][node_num.+included_he] .= true
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

