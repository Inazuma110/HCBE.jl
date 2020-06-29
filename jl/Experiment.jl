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

function experiment(h, hname)
  pyplot()
  # arr = []
  # for b in 0:0.5:1.0
  #   params = Dict("k1"=>2.0, "b"=>b)
  #   k1 = params["k1"]
  #   b = params["b"]
  #   push!(arr, plot_incidence(h, "okapi k1=$k1 b=$b", okapi, params))
  #   params = Dict("k1"=>1.2, "b"=>b)
  #   k1 = params["k1"]
  #   b = params["b"]
  #   push!(arr, plot_incidence(h, "okapi k1=$k1 b=$b", okapi, params))
  # end
  # Plots.plot(arr..., size=(1000, 500))
  # Plots.savefig("./images/$hname-okapi_incidence_matrix.eps")

  # p_arr = []
  # push!(p_arr, plot_incidence(h, "TF-IDF", tfidf))
  # push!(p_arr, plot_incidence(h, "Okapi k1=2.0 b1=1.0", okapi, Dict("k1"=>2.0, "b"=>1.0)))
  # push!(p_arr, plot_incidence(h, "TF", tf))
  # push!(p_arr, plot_incidence(h, "IDF", idf))
  # push!(p_arr, plot_incidence(h, "Random", random_weight))
  #
  # Plots.plot(p_arr..., size=(1000, 500))
  # Plots.savefig("./images/$hname-incidence_matrix.eps")
  #
  # scores = []
  tr_data = Set.([1:100, 101:200, 201:300, 301:400, 401:500])
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, tfidf)
  # push!(scores,scoring.(part_hist, Ref(tr_data), Ref(f1_score)))
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, okapi, Dict("k1"=>2.0, "b"=>1.0))
  # push!(scores,scoring.(part_hist, Ref(tr_data), Ref(f1_score)))
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, tf)
  # push!(scores,scoring.(part_hist, Ref(tr_data), Ref(f1_score)))
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, idf)
  # push!(scores,scoring.(part_hist, Ref(tr_data), Ref(f1_score)))
  # uf,ms, part, part_hist = clustering3(h, 1, modularity, random_weight)
  # push!(scores,scoring.(part_hist, Ref(tr_data), Ref(f1_score)))
  #
  # Plots.plot(scores, label=["TF-IDF" "Okapi k1=2.0 b=1.0" "TF" "IDF" "Random"], lw=5, xlabel="#Added edges", ylabel="F1 measure")
  # Plots.savefig("./images/$hname-f1masures.eps")

  scores2 = []
  for b in 0:0.5:1
    params = Dict("k1"=>1.2, "b"=>b)
    uf,ms, part, part_hist = clustering3(h, 1, modularity, okapi, params)
    push!(scores2,scoring.(part_hist, Ref(tr_data), Ref(f1_score)))
    params = Dict("k1"=>2.0, "b"=>b)
    uf,ms, part, part_hist = clustering3(h, 1, modularity, okapi, params)
    push!(scores2,scoring.(part_hist, Ref(tr_data), Ref(f1_score)))
  end

  Plots.plot(scores2, label=["k1=1.2, b=0" "k1=2.0, b=0" "k1=1.2, b=0.5" "k1=2.0, b=0.5" "k1=1.2, b=1.0" "k1=2.0, b=1.0"], lw=5, xlabel="#Added edges", ylabel="f1 measure")
  Plots.savefig("./images/$hname-okapi_f1masures.eps")

  # cs = []
  # plot_names = []
  # indicator_names = ["Okapi k1=2.0 b=1.0" "TF" "IDF" "Random"]
  # indicators = [okapi, tf, idf, random_weight]
  # for i in 1:5
  #   for j in i+1:5
  #     if i == j continue end
  #     name1 = indicator_names[i]
  #     name2 = indicator_names[j]
  #     ind1 = indicators[i]
  #     ind2 = indicators[j]
  #     push!(cs, curve(h, ind1, ind2))
  #     #         plot_names[p] = "$name1 $name2"
  #     #         p += 1
  #     push!(plot_names, "$name1 $name2")
  #   end
  # end
  #
  # for (name, ind) in zip(indicator_names, indicators)
  #   push!(cs, curve(h, tfidf, ind, f1_score, Dict(), Dict("k1"=>2.0, "b"=>1.0)))
  #   push!(plot_names, "TF-IDF VS $name")
  # end
  #
  # Plots.plot(cs, label=permutedims(plot_names), xlabel="#Added edges", ylabel="f1 measure", xscale=:log10, lw=5)
  # # permutedims(plot_names)
  # Plots.savefig("./images/$hname-f1_curve.eps")
  #
end

function experiment_real(h::Hypergraph, hname)
  ms_hists = []
  ms_hists2 = []
  for i in 1:Int(floor(nhv(h)/100)):nhv(h)
    println(i)
    uf,ms, part, part_hist = clustering3(h, i, my_mod, tf)
    m = my_mod(h, part)
    push!(ms_hists2, m)
  end
  # uf,ms, part, part_hist = clustering3(h, 30000, modularity, tfidf)
  # println(length(part))
  # m = my_mod(h, Set.(1:nhv(h)))

  @save "./tf_ms.jld2" ms_hists2

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

# experiment(rg1, "rg1")
# experiment(rg2, "rg2")
youtube = build_youtube()
experiment_real(youtube, "youtube")
