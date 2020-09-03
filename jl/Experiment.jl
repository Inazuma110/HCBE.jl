using Pkg
Pkg.activate(".")
include("./HypergraphClustering.jl")
include("./IO.jl")
include("./CreateHypergraph.jl")

using SimpleHypergraphs, Random, ProgressBars, Plots, ProgressMeter, JLD2, FileIO, PyPlot, ScikitLearn, GraphPlot
@sk_import metrics : f1_score
@sk_import metrics : accuracy_score
include("HypergraphClustering.jl")
include("./IO.jl")
rg1 = txt2h("./rg.txt")
rg2 = txt2h("./rg2.txt")
rg3 = txt2h("./rg3.txt")
rg4 = txt2h("./rg4.txt")
youtube = build_youtube()


function experiment(h, hname)
  pyplot()
  tr_data = Set.([1:100, 101:200, 201:300, 301:400, 401:500])
  path = "./data/kbs2020/"

  f = open("$path$hname-data", "w")

  ms, ph, bp, ufh = clustering3(h)
  Plots.plot(ms, xlabel="#Added edges in bipartite graph", ylabel="Modularity")
  Plots.savefig("$path$hname-hard.eps")

  cfm = CFModularityCNMLike(1000)
  a, b, hist = findcommunities(h, cfm)
  Plots.plot(hist, xlabel="#Simulations", ylabel="Modularity")
  Plots.savefig("$path$hname-existing.eps")

  ms, eph, bcn, uf = he_clustering(h, 5)
  prob_dists = calc_entropy(h, uf, eph[end])
  for (k, v) in prob_dists
    if !(1.0 in v) println(f, k, ' ', v) end
  end
  avg_ents = avg_entropy.(values(prob_dists))
  arr = []
  for (i, (k, v)) in enumerate(prob_dists)
    push!(arr, (k, avg_ents[i]))
  end
  sort!(arr, rev=true, by=i->i[2])
  println(f, arr)

  Plots.histogram(avg_ents, xlabel="Entropy", ylabel="#Node", bins=10, yscale=:log10)
  Plots.savefig("$path$hname-ent_hist.eps")
  close(f)
  return hist
end

function experiment_real(h::Hypergraph, hname)
  # ms_hists = []
  # fs = [tfidf, okapi, tf, idf, random_weight]
  # tr_data = Set.([1:100, 101:200, 201:300, 301:400, 401:500])
  # for ind in fs
  #   ms_hists2 = []
  #   # for i in 1:Int(floor(nhv(h)/100)):nhv(h)
  #     ms, ph, bp, ufh, scores = clustering3(h, 1, modularity, ind)
  #     score = scoring.(Ref(tr_data), ph, f1_score)
  #     # push!(ms_hists2, score)
  #   # end
  #   # push!(ms_hists, ms_hists2)
  #   push!(ms_hists, score)
  # end
# include("./ClusteringUtil.jl")
# params = Dict("k1"=> 2.0, "b"=>0.75)
  arr = []
  for b in 0:0.5:1.0
    params = Dict("k1"=>1.2, "b"=>b)
    ms, ph, bp, ufh, scores = clustering3(h, 1, my_mod, okapi, params, freq=100)
    push!(arr, ms)
    params = Dict("k1"=>2.0, "b"=>b)
    ms, ph, bp, ufh, scores = clustering3(h, 1, my_mod, okapi, params, freq=100)
    push!(arr, ms)
  end

  @save "./youtube_okapis.jld2" arr

  # push!(ms_hists, ms)
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, okapi)
  # push!(ms_hists, ms)
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, tf)
  # push!(ms_hists, ms)
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, idf)
  # push!(ms_hists, ms)
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, random_weight)
  # push!(ms_hists, ms)
  #
  # Plots.plot(scores, label=["TF-IDF" "Okapi k1=2.0 b=1.0" "TF" "IDF" "Random"], lw=5, xlabel="#added edges", ylabel="f1 measure")
  # Plots.savefig("./images/$hname-f1masures.eps")
end


experiment_real(youtube, "youtube")
# arr = []
# push!(arr, experiment(rg1, "rg1"))
# push!(arr, experiment(rg2, "rg2"))
# push!(arr, experiment(rg3, "rg3"))
# push!(arr, experiment(rg4, "rg4"))
# Plots.plot(arr, labels=["rg1" "rg2" "rg3" "rg4"], xlabel="#Added edges in bipartite graph", ylabel="Modularity", lw=5)
# Plots.savefig("./data/kbs2020/existing_all.eps")
#
