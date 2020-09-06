using Test, HCBE, SimpleHypergraphs

h = Hypergraph{Int64}(6,9)
h[1:3,1] .= 1
h[1:2,2] .= 1
h[2:3,3] .= 1
h[1, 4] = 1
h[3, 4] = 1
h[3:4, 5] .= 1
h[4:6,6] .= 1
h[4:5,7] .= 1
h[5:6, 8] .= 1
h[4, 9] = 1
h[6, 9] = 1
# h[3, 5:9] .= 1

@testset "Clustering algorithm" begin
  @testset "h_HCBE method" begin
    ms, ph, bp, ufh, scores = HCBE.h_HCBE(h, n_cluster=1)
    @test ms == zeros(length(ms))
    @test ph[end] == [Set{Int64}(1:nhv(h))]
    @test length(ph[end]) == 1
    ms, ph, bp, ufh, scores = HCBE.h_HCBE(h, n_cluster=3)
    @test length(ph[end]) == 3
    ms, ph, bp, ufh, scores = HCBE.h_HCBE(h, n_cluster=3, freq=1)
    @test ms != zeros(length(ms))
  end

  @testset "s_HCBE method" begin
    ms, eph, bcn, uf = HCBE.s_HCBE(h, n_cluster=1)
    @test eph[end] == [Set{Int64}(nhv(h)+1:nhe(h)+nhv(h))]
    @test length(eph[end]) == 1
    ms, eph, bcn, uf = HCBE.s_HCBE(h, n_cluster=3)
    @test length(eph[end]) == 3
  end
end

@testset "Clustering util" begin
end
