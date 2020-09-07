"""
  edge

Edge of bipartite graph with star expansion of a hypergraph.

**Arguments**

* from : Vertex number that the edge connects to.  It is a number representing a vertex of the hypergraph.
* to : Vertex number that the edge connects to. It is a number representing a hyperedge of the hypergraph.
* weight : Edge weight.
* id : The number of edges that this structure represents.
"""
mutable struct edge
  from::Int
  to::Int
  weight::Float64
  id::Int
end

"""
  edge_comp(a::edge, b::edge)

Edge comparison function.
Sort in descending order of weight.
If the weights are the same, they are compared according to the IDs of the edges.
"""
function edge_comp(a::edge, b::edge)
  if a.weight == b.weight
    return a.id < b.id
  end
  return a.weight < b.weight
end

"""
partition(uf::UnionFind{Int},
          to::Int=length(uf.parent),
          from::Int=1,
          )::Vector{Set{Int}}

Compute the clustering result from the disjoint set.

**Arguments**

* uf : The disjoint set obtained by clustering.
* to : Disjoint sets up to `to` are used. For example, if you want to get a vertex-only clustering result,
specify the number of vertices.
* from : Use a disjoint set from `from`.
"""
function partition(uf::UnionFind{Int64},
                   to::Int=length(uf.parent),
                   from::Int=1
                  )::Vector{Set{Int}}
  d = Dict()
  for (i, elem) in (enumerate(uf.parent[from:to]))
    cluster_num = root(uf, from-1+i)
    d[cluster_num] = get(d, cluster_num, Set())
    push!(d[cluster_num], from-1+i)
  end
  d = [(length(v), Set{Int64}(v)) for (k, v) in d]
  sort!(d, by=x->x[1], rev=true)
  d = Vector([s for (l, s) in d])
  return d
end

"""
  my_mod(h::Hypergraph, part::Array{Set{any}})

Computes the modularity when `h` is divided by `part`.
"""
function my_mod(h::Hypergraph, part)
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

  (ml - mr) / nhe(h)
end

"""
  build_bg(h::Hypergraph, weighted_f::Function=tfidf, param=Dict())

Star expand `h` and construct a bipartite graph. Each edge is weighted by `weighted_f`.
"""
function build_bg(h::Hypergraph, weighted_f=tfidf, params=Dict())
  edges = Set()
  params["dl_ave"] = 0

  for he in 1:nhe(h)
    params["dl_ave"] += length(getvertices(h, he))
  end
  params["dl_ave"] /= nhe(h)

  id = 1
  @showprogress 1 "computing..." for node in (1:nhv(h))
    hes = gethyperedges(h, node)
    hes = [key for (key, val) in hes]
    for he in hes
      w = weighted_f(h, he, node, params)
      # edgeはnode, edge, tfidfの構造体
      push!(edges, edge(node, nhv(h)+he, w, id))
      id += 1
    end
  end
  edges = [i for i in edges]
  sort!(edges, rev=true, lt=edge_comp)
  return edges
end

function part2is_samecluster(part)
  node_num = sum(length.(part))
  samecluster = Dict([i => Dict([j => 0 for j in i+1:node_num]) for i in 1:node_num])

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


function plot_incidence(h::Hypergraph, name="", weighted_f=tfidf, params=Dict())
  # pyplot()
  # gr()
  bg = build_bg(h, weighted_f, params)
  plot_arr::Array{Tuple{Int64, Int64}} = [(i.to-nhv(h), i.from) for i in bg]
  maxw = maximum([i.weight for i in bg])
  minw = minimum([i.weight for i in bg])
  f(w) = sqrt((w-minw) / (maxw-minw))
  color = [RGB(f(i.weight), f(i.weight), f(i.weight)) for i in bg]
  # color = [RGB(0, f(i.weight), 1.) for i in bg]
  return Plots.scatter(
                plot_arr,
                markersize=6,
                markerstrokewidth=0,
                color=color,
                xlabel="Hyperedges",
                ylabel="Node",
                label="",
                title=name,
               )
end

function plot_incidence2(h::Hypergraph, arr)
  # pyplot()
  # gr()
  bg = arr
  # bg = build_bg(h, weighted_f, params)
  plot_arr::Array{Tuple{Int64, Int64}} = [(i.to-nhv(h), i.from) for i in bg]
  maxw = length(bg)
  minw = 1
  f(w) = sqrt((w-minw) / (maxw-minw))
  color = [RGB(f(i), f(i), f(i)) for (i, e) in enumerate(bg)]
  # println(color)
  # color = [RGB(0, f(i.weight), 1.) for i in bg]
  return Plots.scatter(
                plot_arr,
                markersize=6,
                markerstrokewidth=0,
                color=color,
                xlabel="hyperedges",
                ylabel="node",
                label="",
               )
end

function scoring(tr_part, pred_part, scoring_f)
  act_sc = part2is_samecluster(tr_part)
  pred_sc = part2is_samecluster(pred_part)
  return scoring_f(act_sc, pred_sc)
end

function curve(h::Hypergraph, f1, f2, similar_f=f1_score, params1=Dict(), params2=Dict())
  bg1 = build_bg(h, f1, params1)
  bg2 = build_bg(h, f2, params2)
  # bg1 = map(i->i.id, bg1[1:Int(floor(length(bg1)/100)):end])
  # bg2 = map(i->i.id, bg2[1:Int(floor(length(bg1)/100)):end])
  bg1 = map(i->i.id, bg1)
  bg2 = map(i->i.id, bg2)
  n = length(bg1)

  # arr1 = [0 for i in 1:n]
  # println(length(arr1[bg1]))
  res = []
  @showprogress 1 "computing..."  for i in 1:n
    arr1 = [0 for i in 1:n]
    arr2 = [0 for i in 1:n]
    arr1[bg1[1:i]] .= 1
    arr2[bg2[1:i]] .= 1
    # arr2 = bg2[1:i]
    push!(res, similar_f(arr1, arr2))
  end

  return res
end

function disp_basic(h::Hypergraph)
  node_degree_dist = [0 for i in 1:nhe(h)]
  he_degree_dist = [0 for i in 1:nhv(h)]
  for node in 1:nhv(h)
    size = length(gethyperedges(h, node))
    node_degree_dist[size] += 1
  end
  for he in 1:nhe(h)
    size = length(getvertices(h, he))
    he_degree_dist[size] += 1
  end
  # println(node_degree_dist)
  # println(he_degree_dist)
  p = Plots.bar(node_degree_dist, show=true)
  display(p)
  p = Plots.bar(he_degree_dist, show=true)
  display(p)
end

function epart2cluster(h::Hypergraph, epart)
  part::Array{Set{Int64}} = []
  @showprogress for ecluster in epart
    p = Set([])
    for he in ecluster
      p = union(p, keys(getvertices(h, he-nhv(h))))
    end
    push!(part, p)
  end
  return part
end

function calc_entropy(h::Hypergraph, uf, p; weighted_f=tfidf, params=Dict())
  cluster_set = Set()
  for i in nhv(h)+1:nhv(h)+nhe(h) push!(cluster_set, root(uf, i)) end
  cluster_nums = Dict([val => i for (i, val) in enumerate(cluster_set)])
  entropies = Dict([i => [.0 for j in 1:length(cluster_nums)] for i in 1:nhv(h)])
  bg = build_bg(h, weighted_f)
  @showprogress 1 "computing.." for e in bg
    cluster_num = cluster_nums[root(uf, e.to)]
    entropies[e.from][cluster_num] += e.weight
  end
  for e in bg
    entropies[e.from] ./= sum(entropies[e.from])
  end
  return entropies
end

function avg_entropy(prob_dist)
  avg = 0
  for prob in prob_dist
    if prob == 0 continue end
    avg += prob * log(length(prob_dist), prob)
  end
  avg *= -1
  if avg == -0.0 return 0
  else return avg end
end

function prob_dist2part(prob_dists, threshold=0.5)
  part::Array{Set{Int64}} = [Set() for i in 1:length(prob_dists[1])]
  @showprogress 1 "computing..." for (node, prob_dist) in prob_dists
    for (i, prob) in enumerate(prob_dist)
      if prob > threshold push!(part[i], node) end
    end
  end
  return part
end

function cluster_visualization(h::Hypergraph, p)
  colors = ["red" "blue" "yellow" "green" "brown" "cyan" "magenta" "pink" "white" "gray" ]
  ec = ["gray" for i in nhe(h)]
  c = ["" for i in 1:nhv(h)]
  sort!(p, rev=true, by=i->length(i))

  for (i, cluster) in enumerate(p)
    for node in cluster
      c[node] = length(cluster) == 1 ? "black" : colors[i]
    end
  end

  SimpleHypergraphs.draw(h,
                         HyperNetX,
                         nodes_kwargs=Dict(["facecolors"=>c]),
                         edges_kwargs=Dict(["edgecolors"=>ec]),
                         layout_kwargs=Dict(["seed"=>0]),
                         # with_node_labels=false
                        )

  return c
end

function disp_cluster_bias(p)
  Plots.pie(length.(p), label="")
end

function h2correlation(h::Hypergraph, f1, f2)
  okapi_e = build_bg(h, f1)
  tfidf_e = build_bg(h, f2)

  rank1 = Dict([i => Dict([j => 0 for j in 1:nhe(h)]) for i in 1:nhv(h)])
  rank2 = Dict([i => Dict([j => 0 for j in 1:nhe(h)]) for i in 1:nhv(h)])
  for (i, e) in enumerate(okapi_e)
    rank1[e.from][e.to-nhv(h)] = i
  end

  for (i, e) in enumerate(tfidf_e)
    rank2[e.from][e.to-nhv(h)] = i
  end

  arr1 = Array{Int64}([])
  arr2 = Array{Int64}([])
  @showprogress 1 "computing..." for node in 1:nhv(h)
    for he in 1:nhe(h)
      if rank1[node][he] != 0 push!(arr1, rank1[node][he]::Int64) end
      if rank2[node][he] != 0 push!(arr2, rank2[node][he]::Int64) end
    end
  end

  return (arr1, arr2)
end

# Or later, feature functoins
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
# function testmod(h::Hypergraph, method::CFModularityCNMLike, fx)
#     ha = HypergraphAggs(h)
#     best_modularity = 0
#     comms = [Set(i) for i in 1:nhv(h)]
#     mod_history = Vector{Float64}(undef, method.reps)
#     for rep in 1:method.reps
#         he = rand(1:nhe(h))
#         vers = collect(keys(getvertices(h, he)))
#         if length(vers) == 0
#             continue
#         end;
#         c = deepcopy(comms)
#         i0 = fx(c, vers)
#         max_i = length(c)
#         i_cur = i0
#         while i_cur < max_i
#             i_cur += 1
#             if length(intersect(c[i_cur],vers)) > 0
#                 union!(c[i0], c[i_cur])
#                 c[i_cur]=c[max_i]
#                 max_i += -1
#             end
#                 println(c)
#         end
#         resize!(c,max_i)
#         m = my_mod(h, c)
#         if m > best_modularity
#             best_modularity = m
#             comms = c
#         end
#         mod_history[rep] = best_modularity
#     end
#     return (bm=best_modularity, bp=comms, mod_history=mod_history)
# end
