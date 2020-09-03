mutable struct UnionFind{T <: Integer}
  parent:: Vector{T}  # parent[root] is the negative of the size

  function UnionFind{T}(nodes::T) where T<:Integer
    if nodes <= 0
      throw(ArgumentError("invalid argument for nodes: $nodes"))
    end

    parent = -ones(T, nodes)
    new{T}(parent)
  end
end

UnionFind(nodes::Integer) = UnionFind{typeof(nodes)}(nodes)

function root(uf::UnionFind{T}, x::T)::T where T<:Integer
  if uf.parent[x] < 0
    return x
  else
    # uf.parent[x] = root{T}(uf, uf.parent[x])
    # return uf.parent[x]
    return uf.parent[x] = root(uf, uf.parent[x])
  end
end

function issame(uf::UnionFind{T}, x::T, y::T)::Bool where T<:Integer
  return root(uf, x) == root(uf, y)
end

function size(uf::UnionFind{T}, x::T)::T where T<:Integer
  return -uf.parent[root(uf, x)]
end

function unite!(uf::UnionFind{T}, x::T, y::T)::Bool where T<:Integer
  x = root(uf, x)
  y = root(uf, y)
  if x == y
    return false
  end
  if uf.parent[x] > uf.parent[y]
    x, y = y, x
  end
  # unite smaller tree(y) to bigger one(x)
  uf.parent[x] += uf.parent[y]
  uf.parent[y] = x
  return true
end
