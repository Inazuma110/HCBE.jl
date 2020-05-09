using Pkg
using Random, SimpleHypergraphs

# clusters cluster num
# npcs node per cluster 1つあたりの頂点数
function create_hypergraph(npcs, hepcs, noise_num, he_rate=0.5, noise_rate=he_rate/length(npcs))
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

  for i in 1:noise_num
    he = randsubseq(1:node_num, noise_rate)
    he = Dict([i => 1 for i in he])
    add_hyperedge!(h, vertices=he)
  end

  h, training_data
end

function visu(cre)
  sca::Array{Tuple{Int64, Int64}} = []
  for i in 1:nhe(cre)
    v = [k for (k,v) in getvertices(cre, i)]
    for j in v
      push!(sca, (i, j))
    end
  end
  sca
end
