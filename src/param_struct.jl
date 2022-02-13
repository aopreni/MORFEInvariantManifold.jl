"""
Overview: parametrisation structure management
"""


"""
> mutable struct Parametrisation
- nc : number of monomial combinations
- W : mappings [Υ, Ψ]
- f : reduced dynamics
- nlr : nonlinear right hand side
- comb : combinations the system is solved for
- conj : conjugacy relations
- cmap : mapping index combinations
"""
mutable struct Parametrisation
  nc::Int64
  W::Matrix{ComplexF64}
  f::Matrix{ComplexF64}
  nlr::Matrix{ComplexF64}
  comb::Matrix{Int64}
  conj::Array{Int64}
  cmap::Array{Int64}
  #
  Parametrisation() = new()
  #
end


"""
> initialize_parametrisation!(param,order,ndofs,neq)
It initializes a generic order of the asymptotic expansion
- param : parametrisation
- order : order
- ndofs : number of degrees of freedom of the reduced model
- neq   : number of equations
"""
function initialize_parametrisation!(param::Parametrisation,p::Int64,
                                     ndofs::Int64,neq::Int64,veps=0)
  #
  nperm = ndofs^p
  nc  = floor(Int64,compute_ncomb(p,ndofs))
  param.nc = nc
  if (p>0)
    comb = zeros(Int64,p,nc)
  else
    comb = zeros(Int64,1,nc)
  end
  conj = zeros(Int64,nc)
  cmap = zeros(Int64,nperm)
  #
  if (p>1)
    #
    ith_c = 1
    combs = zeros(Int64,p)
    sorted_combs = zeros(Int64,p)
    #
    ith_c = fill_comb!(comb,cmap,ndofs,p,combs,ith_c,1)
    generate_pmap!(comb,cmap,ndofs,p,combs,sorted_combs,1)
    generate_conj!(conj,comb,cmap,ndofs,p,combs,nc)
    #
  elseif (p==1)
    nm = Int(ndofs/2)
    for i = 1:nm
      conj[i] = i
      conj[i+nm] = -i
    end
    for i = 1:ndofs
      comb[p,i] = i
      cmap[i] = i
    end
  else
    conj[1] = 1
    comb[1,1] = 0
    cmap[1] = 1
  end
  #
  if (veps==1 && p>0)
    [conj[i]=i for i in 1:p]
  end
  #
  param.comb = comb
  param.conj = conj
  param.cmap = cmap
  #
  param.W   = zeros(ComplexF64,neq*2,nc)
  param.nlr = zeros(ComplexF64,neq,nc)
  param.f   = zeros(ComplexF64,ndofs,nc)
  #
  return nothing
end


"""
> generate_conj!(conj,comb,cmap,ndofs,p,combs,nc)
It maps entries to their conjugate
- conj : conjugate
- comb : matrix that stores the powers of the monomials
- cmap : total combinations map
- ndofs : number of degrees of freedom
- p : order of the asymptotic expansion
- combs : combinations
- nc : number of combinations
"""
function generate_conj!(conj,comb,cmap,ndofs,p,combs,nc)
  t = 1
  nm = Int(ndofs/2)
  for i = 1:nc
    #
    if (conj[i]!=0)
      continue
    end
    #
    for j = 1:p
      combs[j] = comb[j,i]
    end
    #
    conj[i] = t
    #
    for j = 1:p
      if (combs[j]>nm)
        combs[j] -= nm
      else
        combs[j] += nm
      end
    end
    sort!(combs)
    #
    pos1 = find_index_global(combs,ndofs,p)
    #
    cm = cmap[pos1]
    if (cm>0 && conj[cm]==0)
      conj[cm] = -i
    end
    t += 1
  end
  #
  return nothing
end


"""
> generate_pmap!(comb::Matrix{Int64},cmap::Array{Int64},ndofs::Int64,p::Int64,combs::Array{Int64},sorted_combs::Array{Int64},counter::Int64)
It computes combinations maps to account for indexes permutations.
- comb : matrix that stores the powers of the monomials
- cmap : total combinations map
- ndofs : number of degrees of freedom
- p : order of the asymptotic expansion
- combs : combinations
- sorted_combs : sorted combinations
- counter : recursion depth
"""
function generate_pmap!(comb::Matrix{Int64},cmap::Array{Int64},ndofs::Int64,p::Int64,
                        combs::Array{Int64},sorted_combs::Array{Int64},counter::Int64)
  for i = 1:ndofs
    combs[counter] = i
    if (counter<p)
      counter += 1
      generate_pmap!(comb,cmap,ndofs,p,combs,sorted_combs,counter)
      counter -= 1
    else
      pos1 = find_index_global(combs,ndofs,p)
      sorted_combs = sort(combs)
      pos2 = find_index_global(sorted_combs,ndofs,p)
      if !(pos1==pos2)
        cmap[pos1] = -cmap[pos2]
      end
    end
  end
end


"""
> fill_comb!(comb::Matrix{Int64},cmap::Array{Int64},ndofs::Int64,p::Int64,combs::Array{Int64},ith_c::Int64,counter::Int64)
It fills the comb matrix with the associated exponentials
- comb : matrix that stores the powers of the monomials
- cmap : total combinations map
- ndofs : number of degrees of freedom
- p : order of the asymptotic expansion
- combs : combinations
- ith_c : ith-combination
- counter : recursion depth
"""
function fill_comb!(comb::Matrix{Int64},cmap::Array{Int64},ndofs::Int64,
                    p::Int64,combs::Array{Int64},ith_c::Int64,counter::Int64)
  #
  if (counter==1)
    for i = 1:ndofs
      combs[counter] = i
      counter += 1
      ith_c = fill_comb!(comb,cmap,ndofs,p,combs,ith_c,counter)
      counter -= 1
    end
  else
    for i = combs[counter-1]:ndofs
      combs[counter] = i
      if (counter<p)
        counter += 1
        ith_c = fill_comb!(comb,cmap,ndofs,p,combs,ith_c,counter)
        counter -= 1
      else
        for j = 1:p
          comb[j,ith_c] = combs[j]
        end
        pos = find_index_global(combs,ndofs,p)
        cmap[pos] = ith_c
        ith_c += 1
      end
    end
  end
  return ith_c
end



"""
> compute_ncomb(param,order,ndofs,neq)
It computes the monomials combinations that need to be computed
- p : order of the expansion
- ndofs : number of degrees of freedom of the reduced model
"""
function compute_ncomb(p::Int64,ndofs::Int64)
  #
  tmp = p+ndofs-1
  nc = 1
  for i = 1:tmp
    nc = nc*i
    if (i<=p)
      nc /= i
    end
  end
  for i = 1:ndofs-1
    nc /= i
  end
  return nc
  #
end


function free_parametrisation!(param::Parametrisation)
  #
  param.W = zeros(ComplexF64,(1,1))
  param.f = zeros(ComplexF64,(1,1))
  param.nlr = zeros(ComplexF64,(1,1))
  param.comb = zeros(Int64,(1,1))
  param.conj = zeros(Int64,1)
  param.cmap = zeros(Int64,1)
  #
  GC.gc()
  #
  return nothing

end



