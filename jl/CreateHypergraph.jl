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

