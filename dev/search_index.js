var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = HCBE","category":"page"},{"location":"#HCBE","page":"Home","title":"HCBE","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [HCBE]","category":"page"},{"location":"#HCBE.UnionFind","page":"Home","title":"HCBE.UnionFind","text":"UnionFind{T <: Integer}\n\nDeprecated!!   In the future, using DataStructures.jl's DisjointSet.\n\n\n\n\n\n","category":"type"},{"location":"#HCBE.edge","page":"Home","title":"HCBE.edge","text":"edge\n\nEdge of bipartite graph with star expansion of a hypergraph.\n\nArguments\n\nfrom : Vertex number that the edge connects to.  It is a number representing a vertex of the hypergraph.\nto : Vertex number that the edge connects to. It is a number representing a hyperedge of the hypergraph.\nweight : Edge weight.\nid : The number of edges that this structure represents.\n\n\n\n\n\n","category":"type"},{"location":"#HCBE.edge_comp-Tuple{HCBE.edge,HCBE.edge}","page":"Home","title":"HCBE.edge_comp","text":"edge_comp(a::edge, b::edge)\n\nEdge comparison function. Sort in descending order of weight. If the weights are the same, they are compared according to the IDs of the edges.\n\n\n\n\n\n","category":"method"},{"location":"#HCBE.epart2cluster-Tuple{SimpleHypergraphs.Hypergraph,Any}","page":"Home","title":"HCBE.epart2cluster","text":"epart2cluster(h::Hypergraph, epart::Array{Set{Int}})\n\n\n\n\n\n","category":"method"},{"location":"#HCBE.h2txt-Tuple{SimpleHypergraphs.Hypergraph,Any}","page":"Home","title":"HCBE.h2txt","text":"h2txt(h::Hypergraph, fname::AbstractString)\n\nTranslate hypergraph to text file. This is in the following form.\n\nN M\n\nN is number of hypergraph vertices. M is number of hypereedges.\n\n\n\n\n\n","category":"method"},{"location":"#HCBE.my_mod-Tuple{SimpleHypergraphs.Hypergraph,Any}","page":"Home","title":"HCBE.my_mod","text":"my_mod(h::Hypergraph, part::Array{Set{any}})\n\nComputes the modularity when h is divided by part.\n\n\n\n\n\n","category":"method"},{"location":"#HCBE.partition","page":"Home","title":"HCBE.partition","text":"partition(uf::UnionFind{Int},           to::Int=length(uf.parent),           from::Int=1,           )::Vector{Set{Int}}\n\nCompute the clustering result from the disjoint set.\n\nArguments\n\nuf : The disjoint set obtained by clustering.\nto : Disjoint sets up to to are used. For example, if you want to get a vertex-only clustering result,\n\nspecify the number of vertices.\n\nfrom : Use a disjoint set from from.\n\n\n\n\n\n","category":"function"},{"location":"#HCBE.s_HCBE-Tuple{SimpleHypergraphs.Hypergraph}","page":"Home","title":"HCBE.s_HCBE","text":"sHCBE(h::Hypergraph,         ncluster::Int=1,         modularityf=modularity,         weightf::Function=tfidf,         freq::Int=Inf)\n\nClustering the vertices of the h using soft-Hypergraph Clustering based on Bipartite Expansion.\n\nArguments\n\nh : Clustered hypergraph.\nncluster : When the number of clusters reaches `ncluster`, the clustering process is terminated.\nmodularity_f : Modularity function.\nweighted_f : Bipartite graph edge weighting function.\nfreq : Once every freq times, the modularity is calculated.\n\n\n\n\n\n","category":"method"},{"location":"#HCBE.star_expansion","page":"Home","title":"HCBE.star_expansion","text":"starexpansion(h::Hypergraph, weightedf::Function=tfidf, param=Dict())\n\nStar expandion h and construct a bipartite graph. Each edge is weighted by weighted_f.\n\n\n\n\n\n","category":"function"}]
}
