include("./Unionfind.jl")

"""
  h_HCBE(h::Hypergraph,
        n_cluster::Int=1,
        modularity_f=modularity,
        weight_f::Function=tfidf,
        freq::Int=Inf)

Clustering the vertices of the `h` using hard-Hypergraph Clustering based on Bipartite Expansion.

**Arguments**

* h : Clustered hypergraph.
* n_cluster : When the number of clusters reaches `n_cluster`, the clustering process is terminated.
* modularity_f : Modularity function.
* weighted_f : Bipartite graph edge weighting function.
* freq : Once every `freq` times, the modularity is calculated.
"""

function h_HCBE(h::Hypergraph;
                n_cluster::Int=1,
                modularity_f=modularity,
                weighted_f=tfidf,
                params=Dict(),
                freq=Inf,
  )
  uf = UnionFind(nhv(h)+nhe(h))
  best_m = 0
  best_part = []
  cluster_num = nhv(h)
  ufh = []
  dl_ave = 0
  for he in 1:nhe(h)
    dl_ave += length(getvertices(h, he))
  end
  dl_ave /= nhe(h)
  params["dl_ave"] = dl_ave

  edges = HCBE.star_expansion(h, weighted_f, params)
  p::Array{Set{Int}} = Set.(1:nhv(h))
  ep = Set.(nhv(h)+1:nhv(h)+nhe(h))
  bcn = -1

  @showprogress 1 "computing..." for (step, edge) in (enumerate(edges))
    node = edge.from
    he = edge.to
    weight = edge.weight
    node_root = root(uf, node)
    he_root = root(uf, he)
    cluster_size = size(uf, node)

    if !issame(uf, node, he)
      unite!(uf, node, he)
      # heのルートがnodeなら
      if he_root <= nhv(h) cluster_num -= 1 end

      if step % freq == 0
        push!(ufh, deepcopy(uf))
      end
    end


    if cluster_num <= n_cluster break end
  end

  if freq != Inf return best_part, ufh
  else return partition(uf, nhv(h)), ufh end
end

"""
  s_HCBE(h::Hypergraph,
        n_cluster::Int=1,
        modularity_f=modularity,
        weight_f::Function=tfidf,
        freq::Int=Inf)

Clustering the vertices of the `h` using soft-Hypergraph Clustering based on Bipartite Expansion.

**Arguments**

* h : Clustered hypergraph.
* n_cluster : When the number of clusters reaches `n_cluster`, the clustering process is terminated.
* modularity_f : Modularity function.
* weighted_f : Bipartite graph edge weighting function.
* freq : Once every `freq` times, the modularity is calculated.
"""
function s_HCBE(h::Hypergraph;
                n_cluster=1,
                modularity_f=modularity,
                weighted_f=tfidf,
                params=Dict(),
                freq=1)
  uf = UnionFind(nhv(h)+nhe(h))
  m = 0
  best_m = -1
  ms = []
  best_part = []
  part_hist = []
  epart_hist = []
  cluster_dict = Dict(1 => nhv(h))
  cluster_num = nhv(h)
  dl_ave = 0
  for he in 1:nhe(h)
    dl_ave += length(getvertices(h, he))
  end
  dl_ave /= nhe(h)
  params["dl_ave"] = dl_ave

  edges = star_expansion(h, weighted_f, params)
  p = Set.(1:nhv(h))
  ep = Set.(nhv(h)+1:nhv(h)+nhe(h))
  bcn = -1
  count = 0


  @showprogress 1 "computing..." for (step, edge) in (enumerate(edges))
    node = edge.from
    he = edge.to
    weight = edge.weight
    node_root = root(uf, node)
    he_root = root(uf, he)
    cluster_size = size(uf, node)

    if !issame(uf, node, he)
      unite!(uf, node, he)
      # heのルートがnodeなら
      if he_root <= nhv(h) cluster_num -= 1 end

      ep = partition(uf, length(uf.parent), nhv(h)+1)
      if step % freq == 0
        # m = modularity(h, epart2cluster(h, ep))
        if m > best_m
          best_m = m
          best_part = p
          bcn = cluster_num
        end
      end
    end

    push!(ms, m)
    push!(epart_hist, ep)

    if length(ep) <= n_cluster break end
  end

  p = epart2cluster(h, ep)
  return ms, p, bcn, uf
end
