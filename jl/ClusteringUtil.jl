using Pkg, SimpleHypergraphs, Random,  ProgressBars, Plots, LightGraphs, ProgressMeter
include("./Unionfind.jl")

mutable struct edge
  from::Int
  to::Int
  weight::Float64
end

function edge_comp(a::edge, b::edge)
  if a.weight == b.weight
    if a.from == b.from return a.to < b.to end
    return a.from < b.from
  end
  return a.weight < b.weight
end

function gini(h::Hypergraph, cluster_dict)
  k = 0
  sum = 0
  for (key, val) in cluster_dict
    k += val
    sum += key * val
  end
  if k == 1 return (1, 1) end

  gini = 0

  for (i, val1) in cluster_dict
    for (j, val2) in cluster_dict
      gini += abs(i - j) * (val1 * val2)
    end
  end
  gini /= 2

  gini = (gini / binomial(k, 2)) / (sum / k)
  gini /= 2

  return (gini, k)
end

function gini2(h::Hypergraph, cluster_dict)
  k = 0
  sum = 0

  for (key, val) in cluster_dict
    if(key == 1) continue end
    k += val
    sum += key * val
  end
  if k <= 1 return (1, 1) end

  gini = 0

  for (i, val1) in cluster_dict
    for (j, val2) in cluster_dict
      if(i == 1 || j == 1) continue end
      gini += abs(i - j) * (val1 * val2)
    end
  end
  gini /= 2

  gini = (gini / binomial(k, 2)) / (sum / k)
  gini /= 2

  return (gini, k)
end

function jaccard(s1, s2)::Float64
  return length(intersect(s1, s2)) / length(union(s1, s2))
end

function simpson(s1, s2)::Float64
  return length(intersect(s1, s2)) / min(length(s1), length(s2))
end

function dice(s1, s2)::Float64
  return (2 * length(intersect(s1, s2))) / (length(s1) + length(s2))
end

function okapi(h::Hypergraph, he::Int, node::Int, dl_ave, k1=2.0, b=0.75)::Float64
  tf = 1.0 / length(getvertices(h, he))
  idf = log(((nhe(h)-length(gethyperedges(h, node)))+0.5) / (length(gethyperedges(h, node))+0.5))
  dl = length(getvertices(h, he))
  v = (idf * tf * (k1+1)) / ((tf * k1) * (1 - b + b * (dl/dl_ave)))
  return v
end

function tfidf(h::Hypergraph, he::Int, node::Int, dl_ave)::Float64
  tf = 1.0 / length(getvertices(h, he))

  # vertices = keys(getvertices(h, he))
  # tf = 1.0 / sum(length.(gethyperedges.(Ref(h), vertices)))
  idf = log(nhe(h) / length(gethyperedges(h, node))) + 1
  return tf*idf
end

function general_weight(h::Hypergraph, he::Int, node::Int, dl_ave)::Float64
  return 1.0 / length(getvertices(h, he))
end

function random_weight(h::Hypergraph=1, he=1, node=1, dl_ave=1)::Float64
  return rand(1)[1]
end


function partition(uf::UnionFind{Int64}, size=length(uf.parent))::Vector{Set{Int}}
  d = Dict()
  for (i, elem) in (enumerate(uf.parent[1:size]))
    cluster_num = root(uf, i)
    d[cluster_num] = get(d, cluster_num, Set())
    push!(d[cluster_num], i)
  end
  d = [(length(v), Set{Int64}(v)) for (k, v) in d]
  sort!(d, by=x->x[1], rev=true)
  d = Vector([s for (l, s) in d])
  # println(d)
  return d
end

function my_mod(h, part)
  ml = 0
  mr = 0
  v = 0
  for i in 1:nhv(h)
    v += length(gethyperedges(h, i))
  end
  dict = Dict([])
  for i in 1:nhe(h)
    dict[length(getvertices(h, i))] = get(dict, length(getvertices(h, i)), 0) + 1
  end

  # vav = 0
  va = Dict([i => 0.0 for i in 1:length(part)])
  for (i, cluster) in enumerate(part)
    for node in cluster
      va[i] += length(gethyperedges(h, node))
      # vav += length(gethyperedges(h, node))
    end
    va[i] /= v
  end

  for (key, val) in dict
    tmp = 0
    for (key2, val2) in va
      tmp += (val2 ^ key)
    end
    tmp *= val
    mr += tmp
  end

  n2c = [-1 for i in 1:nhv(h)]
  for (num, cluster) in enumerate(part)
    for node in cluster
      n2c[node] = num
    end
  end


  for he_i in 1:nhe(h)
    he = [k for (k, v) in getvertices(h, he_i)]
    flag = true
    fnode_cluster = n2c[he[1]]
    for node in he
      if n2c[node] != fnode_cluster
        flag = false
        break
      end
    end
    if flag ml += 1 end
    # if issubset(he, cluster) ml += 1 end
  end

  # println(ml/nhe(h))
  (ml - mr) / nhe(h)
end

function part2is_samecluster(part)
  node_num = sum(length.(part))
  samecluster = Dict([i => Dict([j => 0 for j in i+1:node_num]) for i in 1:node_num])
  # samecluster = [i => Dict(j => 0) for i in 1:node_num for j in 1:node_num]

  for cluster in part
    for n1 in cluster
      for n2 in cluster
        if n2 <= n1 continue end
        samecluster[n1][n2] = true
      end
    end
  end

  samecluster = [samecluster[i][j] for i in 1:node_num for j in i+1:node_num]
  return samecluster
end

function build_bg(h::Hypergraph, weighted_f=tfidf)
  edges = Set()
  dl_ave = 0

  for he in 1:nhe(h)
    dl_ave += length(getvertices(h, he))
  end
  dl_ave /= nhe(h)

  @showprogress 1 "computing..." for node in (1:nhv(h))
    hes = gethyperedges(h, node)
    hes = [key for (key, val) in hes]
    for he in hes
      w = weighted_f(h, he, node, dl_ave)
      # edgeはnode, edge, tfidfの構造体
      push!(edges, edge(node, nhv(h)+he, w))
    end
  end
  edges = [i for i in edges]
  sort!(edges, rev=true, lt=edge_comp)
  return edges
end

# 何番目にnoise nodeが追加されたか
function noise_order(edges, noise_indexes)
  orders = []
  for (i, e) in enumerate(edges)
    node = e.from
    he = e.to
    if noise_indexes[node][he] push!(orders, i) end
  end

  return orders
end


function plot_incidence(h::Hypergraph, name="", weighted_f=tfidf)
  # pyplot()
  # gr()
  bg = build_bg(h, weighted_f)
  plot_arr::Array{Tuple{Int64, Int64}} = [(i.to-nhv(h), i.from) for i in bg]
  maxw = maximum([i.weight for i in bg])
  minw = minimum([i.weight for i in bg])
  f(w) = (w-minw) / (maxw-minw)
  color = [RGB(f(i.weight), f(i.weight), f(i.weight)) for i in bg]
  return Plots.scatter(
                plot_arr,
                markersize=6,
                color=color,
                xlabel="node",
                ylabel="hyperedge",
                label="",
                title=name,
               )
end

function scoring(tr_part, pred_part, scoring_f)
  act_sc = part2is_samecluster(tr_part)
  pred_sc = part2is_samecluster(pred_part)
  return scoring_f(act_sc, pred_sc)
end
