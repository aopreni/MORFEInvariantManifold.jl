"""
Overview: collections of routines for O(n) access to data
"""

"""
> find_index_global(combs,ndofs,p)
It computes the position in the cmap array of a given combination
- combs : total combinations
- ndofs : number of degrees of freedom
- p : order of the expansion
"""
function find_index_global(combs,ndofs,p)
  #
  # It computes the position in the cmap array of a given combination
  #
  pos = 1
  for i = 1:p
    pos += (combs[i]-1)*ndofs^(p-i)
  end
  return pos
end


"""
> find_index_map(j::Int64,k::Int64,s::Int64,p::Int64,ndofs::Int64,index_comb::Array{Int64})
It finds the position of a map for a certain index_comb during time-derivatives residual calculation.
- j : counter
- k : counter
- s : counter
- p : expansion order
- ndofs : number of degrees of freedom
- index_comb : reference monomials exponentials
"""
function find_index_map(j::Int64,k::Int64,s::Int64,p::Int64,
                        ndofs::Int64,index_comb::Array{Int64})
  pos = 1
  @inbounds @simd for i = 1:k
    pos += (index_comb[i]-1)*ndofs^(p-j-i+1)
  end
  @inbounds pos += (s-1)*ndofs^(p-j-k)
  @inbounds @simd for i = k+j+1:p
    pos += (index_comb[i]-1)*ndofs^(p-i)
  end
  return pos
end



"""
> find_index_rdyn(j::Int64,k::Int64,ndofs::Int64,index_comb::Array{Int64})
It finds the position of a reduced dynamics vector for a certain index_comb during time-derivatives residual calculation.
- j : counter
- k : counter
- ndofs : number of degrees of freedom
- index_comb : reference monomials exponentials
"""
function find_index_rdyn(j::Int64,k::Int64,ndofs::Int64,index_comb::Array{Int64})
  pos = 1
  @inbounds @simd for i = k+1:k+j
    pos += (index_comb[i]-1)*ndofs^(k+j-i)
  end
  return pos
end




"""
> find_position_g(j::Int64,p::Int64,ndofs::Int64,index_comb::Array{Int64})
It finds the maps position to compute their contribution to the quadratic nonlinearities.
- j : counter
- p : order
- ndofs : number of degrees of freedom
- index_comb : reference monomials exponentials
"""
function find_position_g(j::Int64,p::Int64,ndofs::Int64,
                          index_comb::Array{Int64})
  pos1 = 1
  pos2 = 1
  @inbounds @simd for i = 1:j
    pos1 += (index_comb[i]-1)*ndofs^(j-i)
  end
  @inbounds @simd for i = j+1:p
    pos2 += (index_comb[i]-1)*ndofs^(p-i)
  end
  return pos1,pos2
end





"""
> find_position_g(j::Int64,k::Int64,p::Int64,ndofs::Int64,index_comb::Array{Int64})
It finds the maps position to compute their contribution to the cubic nonlinearities.
- j : counter
- k : counter
- p : order
- ndofs : number of degrees of freedom
- index_comb : reference monomials exponentials
"""
function find_position_h(j::Int64,k::Int64,p::Int64,
                          ndofs::Int64,index_comb::Array{Int64})
  #
  # It finds the maps position to compute their contribution to the 
  # cubic nonlinearities.
  #
  pos1 = 1
  pos2 = 1
  pos3 = 1
  @inbounds @simd for i = 1:j
    pos1 += (index_comb[i]-1)*ndofs^(j-i)
  end
  @inbounds @simd for i = 1:k
    pos2 += (index_comb[j+i]-1)*ndofs^(k-i)
  end
  @inbounds @simd for i = j+k+1:p
    pos3 += (index_comb[i]-1)*ndofs^(p-i)
  end
  return pos1,pos2,pos3
end
