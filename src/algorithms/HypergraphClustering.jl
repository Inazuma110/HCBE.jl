include("./Unionfind.jl")

function clustering(h::Hypergraph, gini_func::Function=gini, order=shuffle!([i for i in 1:nhe(h)]))
  ginis = []
  res = [-1 for i in 1:nhv(h)]
  cluster_dict = Dict([1 => nhv(h)])
  uf = UnionFind(nhv(h))
  euf = UnionFind(nhe(h))
  done = Set()
  ks = []
  he_visited = [false for i in 1:nhe(h)]
  for (i, he_i) in tqdm(enumerate(order))
    he_visited[he_i] = true
    nodes = getvertices(h, he_i)
    nodes = [node for node in keys(nodes)]
    if isempty(nodes) continue end

    first = nodes[1]

    for node in nodes[2:end]
      first_root = root(uf, first)
      node_root = root(uf, node)
      if !issame(uf, first, node)
        cluster_dict[size(uf, first)] -= 1
        if cluster_dict[size(uf, first)] == 0 delete!(cluster_dict, size(uf, first)) end
        cluster_dict[size(uf, node)] -= 1
        if cluster_dict[size(uf, node)] == 0 delete!(cluster_dict, size(uf, node)) end
        unite!(uf, first, node)
        cluster_dict[size(uf, first)] = get(cluster_dict, size(uf, first), 0) + 1
      end
      if(first_root != root(uf, first) && res[first] == -1) res[first] = i end
      if(node_root != root(uf, node) && res[node] == -1) res[node] = i end
    end

    (giniv, k) = gini_func(h, cluster_dict)
    push!(ginis, giniv)
    push!(done, he_i)
    push!(ks, k)
  end
  return (ginis, ks, uf::UnionFind)
end


function clustering2(h::Hypergraph, indicator::Function=jaccard, order=shuffle!([i for i in 1:nhe(h)]))
  visited = [false for i in 1:nhe(h)]
  euf = UnionFind(nhe(h))
  inds = [-1.0 for i in 1:nhe(h)]
  densities = [(-1.0, i) for i in 1:nhe(h)]
  he_nums = [1 for i in 1:nhe(h)]
  # 各クラスタがどのノードを持っているか
  p_nodes = Dict(Set())
  # ここで各エッジのノードを入れる。
  for edge in tqdm(1:nhe(h))
    p_nodes[edge] = get(p_nodes, edge, Set([node for (node, val) in getvertices(h, edge)]))
    edges = Set()
    nodes = [node for (node, v) in getvertices(h, edge)]
    edge_num = 0
    for node in nodes
      edge_num += length(gethyperedges(h, node))
    end
    # density = edge_num / length(nodes)
    density = edge_num
    density = (isnan(density)) ? 0 : density
    densities[edge] = (density, edge)
  end
  sort!(densities, rev=true)
  order = [edge for (density, edge) in densities]

  # 順に回す
  for (step, he) in tqdm(enumerate(order))
    if visited[he] continue end
    he_root = root(euf, he)
    root_dense = he_nums[he_root] / length(p_nodes[he_root])
    # ランダムに選んだエッジが持つノードを列挙
    nodes = getvertices(h, he)
    nodes = [node for (node, val) in nodes]
    max_ind = inds[he_root]

    # 更にそのノードが持つハイパーエッジを列挙
    edges = Set()
    for node in nodes
      for e in [edge for (edge, val) in gethyperedges(h, node)]
        push!(edges, e)
      end
    end

    # 各ハイパーエッジごとにJaccardなどの指標を計算
    # 降順にソート
    hes_indvs = Vector()
    for next_he in edges
      next_nodes = getvertices(h, next_he)
      next_nodes = [node for (node, val) in next_nodes]

      ind_v = indicator(next_nodes, nodes)
      push!(hes_indvs, (ind_v, next_he))
    end
    sort!(hes_indvs, rev=true)

    # 順にマージ
    for (ind_v, next_he) in hes_indvs
      next_he_root = root(euf, next_he)
      next_ind = indicator(p_nodes[he_root], p_nodes[next_he_root])
      nher_dense = he_nums[next_he_root] / length(p_nodes[next_he_root])
      cluster_he = he_nums[he_root] + he_nums[next_he_root]
      next_dens = cluster_he / (length(p_nodes[he_root]) + length(p_nodes[next_he_root]))

      # 最大指標より次の指標値のほうが高ければマージ
      # あったほうがいいのかない方がいいのか不明
      if visited[next_he] continue end
      # if next_dens < root_dense continue end
      # if next_dens < nher_dense continue end

      if next_ind >= max_ind
        # クラスタが違えば
        if !issame(euf, he_root, next_he_root)
          unite!(euf, he_root, next_he_root)
          next_cluster = union(p_nodes[he_root], p_nodes[next_he_root])
          p_nodes[he_root] = next_cluster
          p_nodes[next_he_root] = next_cluster
          inds[he_root] = next_ind
          inds[next_he_root] = next_ind
          he_nums[he_root] = cluster_he
          he_nums[next_he_root] = cluster_he
          max_ind = next_ind
          visited[he] = true
          visited[next_he] = true
        end
      # あったほうがいいのかない方がいいのか不明
      # else
      #   break
      end
    end
  end

  return euf
end

# Calc all param
function h_HCBE(h::Hypergraph;
                n_cluster=1,
                modularity_f=modularity,
                weighted_f=tfidf,
                params=Dict(),
                freq=Inf)
  uf = UnionFind(nhv(h)+nhe(h))
  m = 0
  best_m = 0
  best_score = -1
  ms = []
  best_part = []
  part_hist = []
  cluster_dict = Dict(1 => nhv(h))
  cluster_num = nhv(h)
  ufh = []
  dl_ave = 0
  for he in 1:nhe(h)
    dl_ave += length(getvertices(h, he))
  end
  dl_ave /= nhe(h)
  params["dl_ave"] = dl_ave

  edges = build_bg(h, weighted_f, params)
  p::Array{Set{Int}} = Set.(1:nhv(h))
  ep = Set.(nhv(h)+1:nhv(h)+nhe(h))
  bcn = -1
  y = Set.([1:100, 101:200, 201:300, 301:400, 401:500])
  scores = []
  score = 0


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

      p = partition(uf, nhv(h))
      if step % freq == 0
        m = modularity_f(h, p)
        # score = scoring(y, p, f1_score)
        if m > best_m
          best_score = score
          best_part = p
          bcn = cluster_num
        end
      end
    end

    push!(ms, m)
    push!(part_hist, p)
    push!(ufh, copy(uf.parent))
    push!(scores, score)

    if cluster_num <= n_cluster break end
  end

  return ms, part_hist, best_part, ufh, scores
end

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

  edges = build_bg(h, weighted_f, params)
  # dendrogram = [(-1, -1, -1) for i in 1:nhv(h)+nhe(h)]
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

      if step % freq == 0
        ep = partition(uf, length(uf.parent), nhv(h)+1)
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

  return ms, epart_hist, bcn, uf
end


# minimum calc
function clustering4(h::Hypergraph, n_cluster=1, modularity_f=modularity, weighted_f=tfidf)
  uf = UnionFind(nhv(h)+nhe(h))
  cluster_num = nhv(h)
  dl_ave = 0
  for he in 1:nhe(h)
    dl_ave += length(gethyperedges(h, he))
  end
  dl_ave /= nhe(h)
  edges = build_bg(h, weighted_f)

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
    end

    if cluster_num == n_cluster break end
  end

  return partition(uf, nhv(h))
end

function d_clustering(h1::Hypergraph, n_cluster=1, modularity_f=modularity, weighted_f=tfidf)
  h = Hypergraph(copy(h1))
  order = []

  cluster_num = nhv(h)
  part_hist = []
  uf = UnionFind(nhv(h)+nhe(h))
  merged = [[false for j in 1:nhe(h)] for i in 1:nhv(h)]
  dl_ave = 0
  for he in 1:nhe(h)
    dl_ave += length(gethyperedges(h, he))
  end
  dl_ave /= nhe(h)
  edges = build_bg(h, weighted_f)
  @showprogress 1 "computing..." for x in 1:length(edges)

    for (step, edge) in (enumerate(edges))
      node = edge.from
      he = edge.to
      weight = edge.weight
      node_root = root(uf, node)
      he_root = root(uf, he)
      cluster_size = size(uf, node)
      if merged[node][he-nhv(h)]
        continue
      end
      is_connects = Dict([i => issame(uf, he, i) for i in nhv(h)+1:nhv(h)+nhe(h)])

      merged[node][he-nhv(h)] = true
      # push!(order, edge)
      if !issame(uf, node, he)
        unite!(uf, node, he)
        # heのルートがnodeなら
        if he_root <= nhv(h) cluster_num -= 1 end

        is_connects2 = Dict([i => issame(uf, he, i) for i in nhv(h)+1:nhv(h)+nhe(h)])
        if is_connects != is_connects2
          new_vertices = Set(keys(getvertices(h, he-nhv(h))))
          merged_hes = Set(he-nhv(h))
          for i in nhv(h)+1:nhv(h)+nhe(h)
            if is_connects[i] == is_connects2[i] continue end
            push!(merged_hes, i-nhv(h))
            push!.(Ref(new_vertices), Set(keys(getvertices(h, i-nhv(h)))))
          end

          merged_hes = [mhe for mhe in merged_hes]
          sort!(merged_hes, rev=true)
          # for mhe in merged_hes
          #   remove_hyperedge!(h, mhe)
          # end
          # he_ind = add_hyperedge!(h)
          # for i in new_vertices h[i, he_ind] = 1 end

          for e in edges
            if !(e.to-nhv(h) in merged_hes) continue end
            new_idf = length(gethyperedges(h, e.from)) - length(intersect(keys(gethyperedges(h, e.from)), merged_hes))+1
            e.weight = log(nhe(h) / new_idf) + 1
            new_tf = length(merged_hes)
            e.weight /= new_tf
          end
          sort!(edges, rev=true, lt=edge_comp)
          break
        end
      end
      push!(part_hist, partition(uf, nhv(h)))
      push!(order, edge)

      if cluster_num <= n_cluster return uf, part_hist, order end
    end
  end
  return uf, part_hist, order
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

