function jaccard(s1, s2)::Float64
  return length(intersect(s1, s2)) / length(union(s1, s2))
end

function simpson(s1, s2)::Float64
  return length(intersect(s1, s2)) / min(length(s1), length(s2))
end

function dice(s1, s2)::Float64
  return (2 * length(intersect(s1, s2))) / (length(s1) + length(s2))
end

function okapi(h::Hypergraph, he::Int, node::Int, params=Dict("k1"=>2.0, "b"=>0.75))::Float64
  tf = 1.0 / length(getvertices(h, he))
  idf = log(((nhe(h)-length(gethyperedges(h, node)))+0.5) / (length(gethyperedges(h, node))+0.5))
  dl = length(getvertices(h, he))
  dl_ave = params["dl_ave"]
  k1 = haskey(params, "k1") ? params["k1"] : 2.0
  b = haskey(params, "b") ? params["b"] : 0.75
  v = (idf * tf * (k1+1)) / ((tf * k1) * (1 - b + b * (dl/dl_ave)))
  return v
end

function tf(h::Hypergraph, he::Int, node::Int, params=Dict())
  return 1.0 / length(getvertices(h, he))
end

function idf(h::Hypergraph, he::Int, node::Int, params=Dict())
  return log(nhe(h) / length(gethyperedges(h, node))) + 1
end

function tfidf(h::Hypergraph, he::Int, node::Int, params=Dict())::Float64
  return tf(h, he, node) * idf(h, he, node)
end

function transpose_tf(h::Hypergraph, he::Int, node::Int, params=Dict())
  # println(length(gethyperedges(h, node)))
  return 1.0 / length(gethyperedges(h, node))
end

function transpose_idf(h::Hypergraph, he::Int, node::Int, params=Dict())
  return log(nhv(h) / length(getvertices(h, he))) + 1
end

function transpose_tfidf(h::Hypergraph, he::Int, node::Int, params=Dict())
  return transpose_tf(h, he, node) * transpose_idf(h, he, node)
end

function tfidf_dfitf(h::Hypergraph, he::Int, node::Int, params=Dict())
  return tfidf(h, he, node) * transpose_tfidf(h, he, node)
end

function general_weight(h::Hypergraph, he::Int, node::Int, dl_ave)::Float64
  return 1.0 / length(getvertices(h, he))
end

function random_weight(h::Hypergraph=1, he=1, node=1, dl_ave=1)::Float64
  return rand(1)[1]
end

