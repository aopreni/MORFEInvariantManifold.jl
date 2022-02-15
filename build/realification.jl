"""
Overview: collection of subroutines for the Parametrisation
"""


"""
> identity_tangency!(Cp::Parametrisation,nm::Int64,neq::Int64,Φₗᵢₛₜ::Array{Int64},ϕ::Matrix{ComplexF64},ω₀::Array{Float64},ζ₀::Array{Float64})
It fills the first order parametrization using eigenmodes and eigenvalues of the mechanical problem.
- Cp : Parametrisation data structure
- nm : number of master modes
- neq : total number of equations
- Φₗᵢₛₜ : master modes list
- ϕ : modes
- ω₀ : eigenfrequencies
- ζ₀ : damping coefficients
- U : displacement field
"""
function identity_tangency!(Cp::Parametrisation,nm::Int64,neq::Int64,
                            Φₗᵢₛₜ::Array{Int64},ϕ::Matrix{Float64},
                            ω₀::Array{Float64},ζ₀::Array{Float64},
                            U::Field)
  for i = 1:nm
    λ₁ = - ζ₀[i]*ω₀[i] + ω₀[i]*sqrt(1.0-ζ₀[i]^2.0)*im
    λ₂ = - ζ₀[i]*ω₀[i] - ω₀[i]*sqrt(1.0-ζ₀[i]^2.0)*im
    Cp.f[i,i]       = λ₁
    Cp.f[i+nm,i+nm] = λ₂
    @inbounds for j = 1:U.neq
      Cp.W[j,       i]    = ϕ[j,Φₗᵢₛₜ[i]]*λ₁
      Cp.W[j+neq,   i]    = ϕ[j,Φₗᵢₛₜ[i]]
      Cp.W[j,    i+nm]    = ϕ[j,Φₗᵢₛₜ[i]]*λ₂
      Cp.W[j+neq,i+nm]    = ϕ[j,Φₗᵢₛₜ[i]]
    end
  end
  #
  return nothing
  #
end


"""
> realification!(Cp,Wr,fr,p,ndofs,neq,rmat)
It realifies mappings and reduced dynamics for an order p of the parametrisation
- Cp : Parametrisation data structure
- Wr : real-valued mapping
- fr : real-valued reduced dynamics
- p : order of the asymptotic development
- ndofs : number of degrees of freedom
- neq : number of equations
- rmat : realification matrix R
"""
function realification!(Cp,Wr,fr,p,ndofs,neq,rmat)
  #
  # It makes mappings and reduced dynamics real-valued
  #
  Cic = zeros(Int64,p)
  Ric = zeros(Int64,p)
  pm_comb = zeros(Int64,p)
  #
  if (p==0)
    @inbounds Wr[:] = Cp.W[:]
    @inbounds fr[:] = rmat*Cp.f[:]
    return nothing
  end
  #
  for i = 1:Cp.nc
    for j = 1:p
      Cic[j] = Cp.comb[j,i]
    end
    recursive_pm!(Cp,Wr,fr,Cic,Ric,pm_comb,p,ndofs,neq,i,1)
  end
  for i = 1:Cp.nc
    @inbounds fr[:,i] = rmat*fr[:,i]
  end
  return nothing
end


"""
> recursive_pm!(Cp,Wr,fr,Cic,Ric,pm_comb,p,ndofs,neq,comb,cc)
It realifies mappings and reduced dynamics of a given index combination comb of order p.
- Cp : parametrisation data structure
- Wr : real-valued mappings
- fr : real-valued reduced dynamics
- Cic : complex-valued parametrisation index combination of a given monomial
- Ric : real-valued parametrisation index combination of a given monomial
- pm_comb : +- combinations
- p : order of the asymptotic development
- ndofs : number of degrees of freedom
- neq : number of equations
- comb : reference index combination of the realified monomial
- cc : recursion depth
"""
function recursive_pm!(Cp,Wr,fr,Cic,Ric,pm_comb,p,ndofs,neq,comb,cc)
  #
  # Conversion of complex combinations to real-valued combinations.
  #
  nm = Int(ndofs/2)
  bin = [0, Int(ndofs/2)]
  for i = 1:2
    pm_comb[cc] = bin[i]
    if (cc<p)
      cc += 1
      recursive_pm!(Cp,Wr,fr,Cic,Ric,pm_comb,p,ndofs,neq,comb,cc)
      cc -= 1
    else
      pm = 1.0+0im
      #
      @inbounds for j = 1:p
        if (Cic[j]<=bin[2])
          Ric[j] = Cic[j]
        else
          Ric[j] = Cic[j] - bin[2]
        end
      end
      #
      @inbounds for j = 1:p
        if (pm_comb[j]!=0) 
          pm *= 1.0/1.0im
          if (Cic[j]>bin[2])
            pm *= -1.0
          end
        end
      end
      #
      Ric = Ric + pm_comb
      pos = find_index_global(Ric,ndofs,p)
      cm = Cp.cmap[pos]
      cma = abs(cm)
      for j = 1:neq*2
        @inbounds Wr[j,cma] += (pm/2.0^p) * Cp.W[j,comb]
      end
      for j = 1:ndofs
        @inbounds fr[j,cma] += (pm/2.0^p) * Cp.f[j,comb]
      end
    end
  end
  return nothing
end


"""
> init_realification_matrix(ndofs::Int64)
It initializes the matrix used to realify the system
- 
"""
function init_realification_matrix(ndofs::Int64)
  rmat = zeros(ComplexF64,(ndofs,ndofs))
  nm = Int(ndofs/2)
  @inbounds for i = 1:nm
    rmat[i,i] = 1.0
    rmat[i,i+nm] = 1.0
    rmat[i+nm,i] = 1.0im
    rmat[i+nm,i+nm] = -1.0im
  end
  return rmat
end