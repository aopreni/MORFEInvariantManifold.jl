"""
Overview: structures and functions used to manage sparse matrices
nodal connectivity.
"""

"""
> struct NodalConnectivity
- nval : number of nodes + 1
- nreg : total size of the nodal connectivity
- register : bookkeeping nodes to nodes connectivity in compact sparse format
- values : stored node to node connectivity
- maxr : FILL
"""
struct NodalConnectivity
  nval::Int64
  nreg::Int64
  register::Vector{Int64}
  values::Vector{Int64}
  maxr::Int64
end


"""
> NodalConnectivity(mesh::Grid)
NodalConnectivity main constructor.
- mesh : Grid
"""
function NodalConnectivity(mesh::Grid)
  #
  itab = zeros(Int64,(ncv_max*nne_max,mesh.nn))
  kconn = zeros(Int64,mesh.nn)
  inspected = zeros(Bool,mesh.nn)
  # first connectivity build
  for d = 1:mesh.nΩ
    iΩ = mesh.Ω[d]
    for se = 1:iΩ.Sen # iterate over supported element types
      nn = iΩ.Senn[se] #number of nodes per element
      skip = iΩ.eskip[se]
      for e = 1:iΩ.ne[se]
        conn = @view iΩ.e2n[skip+1+(e-1)*nn:skip+e*nn] # gather connectivity
        for k1 = 1:nn
          iv1 = conn[k1]
          for k2 = 1+k1:nn
            iv2 = conn[k2]
            kconn[iv1] += 1
            kconn[iv2] += 1
            itab[kconn[iv1],iv1] = iv2
            itab[kconn[iv2],iv2] = iv1
          end
        end
      end
    end
  end
  # refine
  for iv = 1:mesh.nn
    j = 0
    n = kconn[iv]
    for i = 1:n
      inb = itab[i,iv]
      if (!inspected[inb])
        j += 1
        itab[j,iv] = inb
        inspected[inb] = true
      end
    end
    kconn[iv] = j
    for i = 1:j
      inspected[itab[i,iv]] = false
    end
  end
  # allocate variables
  n = 0
  for iv = 1:mesh.nn
    n += kconn[iv]
  end
  register = Vector{Int64}(undef,mesh.nn+1)
  values = zeros(Int64,n)
  # fill 
  i = 1
  register[i] = i
  for iv = 1:mesh.nn
    i += kconn[iv]
    register[iv+1] = i
    iv1 = register[iv]
    iv2 = register[iv+1]-1
    values[iv1:iv2] = itab[1:kconn[iv],iv]
  end
  #
  return NodalConnectivity(mesh.nn+1, n, register, values, maximum(kconn))
end



"""
> NodalConnectivity(n2n::NodalConnectivity, node::Int64)
It returns a View containing which nodes are connector to node
- n2n : NodalConnectivity structure
- node : node of interest
"""
function get_node_connectivity(n2n::NodalConnectivity,node::Int64)
  #
  return @view n2n.values[n2n.register[node]:n2n.register[node+1]-1]
  #
end