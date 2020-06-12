using SimpleHypergraphs
"""
  h2txt(h::Hypergraph, fname::AbstractString)

Translate hypergraph to text file. This is in the following form.
```
N M
```
N is number of hypergraph vertices.
M is number of hypereedges.

"""
function h2txt(h::Hypergraph, fname)
  out = open(fname, "w")
  print(out, nhv(h), ' ', nhe(h), '\n')
  for he in 1:nhe(h)
    print(out, length(getvertices(h, he)), ' ')
    for node in keys(getvertices(h, he))
      print(out, node, ' ')
    end
    print(out, '\n')
  end
  close(out)
end

function txt2h(fname)
  f = open(fname)
  (n, m) = parse.(Int, split(readline(f), ' '))
  h = Hypergraph{Int}(n, m)
  he_i = 1
  while(!eof(f))
    line = split(readline(f), ' ')
    if length(line) == 0 continue end
    vᵢ= line[1]
    he = parse.(Int, line[2:end-1])
    h[he, he_i] .= 1
    he_i += 1
  end
  close(f)

  return h
end
