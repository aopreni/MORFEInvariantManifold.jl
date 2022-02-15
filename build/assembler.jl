"""
Overview: collections of routines and utilities to
perform matrix assembly
"""


"""
> assembly_forcing_vector!(Cp⁺,mesh,U,M,ϕ,κ_modes,κ_list,κ_phase)
It assemblies the non-autonomous forcing vector

\$ \\mathbf{F} = ∑_{i=1}^{N} κᵢMϕᵢ \$

- Cp⁺ : no-autonomous parametrisation data structure
- mesh : mesh data structure
- U : displacement field
- M : mass matrix
- ϕ : eigenmodes
- κ_modes : list of forced modes
- κ_list : forcing amplitude κᵢ
- κ_phase : phase of the applied forcing
"""
function assembly_forcing_vector!(Cp⁺::Parametrisation,mesh::Grid,U::Field,
                                  M,ϕ,κ_modes,κ_list,κ_phase)
  #
  for (i,κ) in enumerate(κ_modes)
    if (κ_phase[i]=='c')
      Cp⁺.nlr -= (1.0/2.0)*κ_list[i]*M*ϕ[:,κ]
    elseif (κ_phase[i]=='s')
      Cp⁺.nlr -= (1.0/2.0/im)*κ_list[i]*M*ϕ[:,κ]
    end
  end
  #
  return nothing                               
  #
end


"""
> assembly_H_nl!(Cp, entry,Ψ₁, Ψ₂, Ψ₃,mesh, U, mult = 1.0)
It assemblies cubic nonlinearities operator

\$ G(Ψ₁,Ψ₂,Ψ₃) = \\frac{1}{6} ∫_{Ω} γ(Ψ₁,Ψ₂):\\mathcal{A}:γ(Ψ₃,w) + γ(Ψ₁,Ψ₃):\\mathcal{A}:γ(Ψ₂,w) + γ(Ψ₃,Ψ₂):\\mathcal{A}:γ(Ψ₁,w) dΩ \$

- Cp : parametrisation data structure
- entry : entrance of the reference array
- Ψ₁ : mapping
- Ψ₂ : mapping
- Ψ₃ : mapping
- mesh : mesh data structure
- U : displacement field
- mult : integral multiplier
"""
function assembly_H_nl!(Cp::Parametrisation, entry::Int64,
                        Ψ₁, Ψ₂, Ψ₃,
                        mesh::Grid, U::Field, mult = 1.0)
  #
  X = zeros(Float64,nne_max*dim)
  dofs = zeros(Int64,nne_max*dim)
  Fₑ = zeros(ComplexF64,nne_max*dim)
  #
  Ψ₁ₑ = zeros(ComplexF64,nne_max*dim)
  Ψ₂ₑ = zeros(ComplexF64,nne_max*dim)  
  Ψ₃ₑ = zeros(ComplexF64,nne_max*dim) 
  #
  N     = zeros(Float64,nne_max)
  ∂N∂a  = zeros(Float64,(nne_max,dim))
  ∂N∂x  = zeros(Float64,(nne_max,dim))
  Jac   = zeros(Float64,(dim,dim))
  Jac⁻¹ = zeros(Float64,(dim,dim))
  #
  ∇U₁ = zeros(ComplexF64,(dim,dim))
  ∇U₂ = zeros(ComplexF64,(dim,dim))
  ∇U₃ = zeros(ComplexF64,(dim,dim))
  #
  e₁₂ = zeros(ComplexF64,(dim,dim)) 
  e₂₃ = zeros(ComplexF64,(dim,dim)) 
  e₁₃ = zeros(ComplexF64,(dim,dim)) 
  #
  eᵛ₁₂ = zeros(ComplexF64,dim*(dim-1)) 
  eᵛ₂₃ = zeros(ComplexF64,dim*(dim-1)) 
  eᵛ₁₃ = zeros(ComplexF64,dim*(dim-1))
  #
  σᵛ₁₂ = zeros(ComplexF64,dim*(dim-1))
  σᵛ₂₃ = zeros(ComplexF64,dim*(dim-1))
  σᵛ₁₃ = zeros(ComplexF64,dim*(dim-1))
  #
  sym∇ⁿˡ₁ = zeros(ComplexF64,(nne_max*dim,dim*(dim-1)))
  sym∇ⁿˡ₂ = zeros(ComplexF64,(nne_max*dim,dim*(dim-1)))
  sym∇ⁿˡ₃ = zeros(ComplexF64,(nne_max*dim,dim*(dim-1)))
  #
  for iΩ ∈ mesh.Ω
    Dᵢⱼₖₗ = iΩ.mat.Dᵢⱼₖₗ
    for set = 1:iΩ.Sen
      etype = iΩ.Set[set]
      qr = select_quadrature_rule(etype)
      nn = iΩ.Senn[set]
      skip = iΩ.eskip[set]
      for e = 1:iΩ.ne[set]
        conn = @view iΩ.e2n[skip+1+(e-1)*nn:skip+e*nn]
        get_coor!(mesh, conn, X, nn)
        dofs!(U,nn,conn,dofs)
        #
        @inbounds for i = 1:nn*dim
          if (dofs[i]>0)
            Ψ₁ₑ[i] = Ψ₁[dofs[i]]
            Ψ₂ₑ[i] = Ψ₂[dofs[i]]
            Ψ₃ₑ[i] = Ψ₃[dofs[i]]
          else
            Ψ₁ₑ[i] = 0.0
            Ψ₂ₑ[i] = 0.0
            Ψ₃ₑ[i] = 0.0
          end
        end
        #
        integrate_H!(Fₑ,X,Dᵢⱼₖₗ,Ψ₁ₑ,Ψ₂ₑ,Ψ₃ₑ,
                    N,∂N∂a,∂N∂x,Jac,Jac⁻¹,
                    ∇U₁,∇U₂,∇U₃,
                    e₁₂,e₂₃,e₁₃,
                    eᵛ₁₂,eᵛ₂₃,eᵛ₁₃,
                    σᵛ₁₂,σᵛ₂₃,σᵛ₁₃,
                    sym∇ⁿˡ₁,sym∇ⁿˡ₂,sym∇ⁿˡ₃,
                    nn,etype,qr)
        #
        for i = 1:nn*dim
          if (dofs[i]>0)
            @inbounds Cp.nlr[dofs[i],entry] += Fₑ[i]*mult
          end
        end
        #
      end
    end
  end
  #
  return nothing
  #
end



"""
> assembly_G_nl!(Cp, entry,Ψ₁, Ψ₂,mesh, U, mult = 1.0)
It assemblies quadratic nonlinearities operator

\$ G(Ψ₁,Ψ₂) = \\frac{1}{2} ∫_{Ω} γ(Ψ₁,Ψ₂):\\mathcal{A}:ε(w) + γ(Ψ₁,w):\\mathcal{A}:ε(Ψ₂) + γ(w,Ψ₂):\\mathcal{A}:ε(Ψ₁) dΩ \$

- Cp : parametrisation data structure
- entry : entrance of the reference array
- Ψ₁ : mapping
- Ψ₂ : mapping
- mesh : mesh data structure
- U : displacement field
- mult : integral multiplier
"""
function assembly_G_nl!(Cp::Parametrisation, entry::Int64,
                        Ψ₁, Ψ₂,
                        mesh::Grid, U::Field, mult = 1.0)
  #
  X = zeros(Float64,nne_max*dim)
  dofs = zeros(Int64,nne_max*dim)
  Fₑ = zeros(ComplexF64,nne_max*dim)
  #
  Ψ₁ₑ = zeros(ComplexF64,nne_max*dim)
  Ψ₂ₑ = zeros(ComplexF64,nne_max*dim)
  #
  N     = zeros(Float64,nne_max)
  ∂N∂a  = zeros(Float64,(nne_max,dim))
  ∂N∂x  = zeros(Float64,(nne_max,dim))
  Jac   = zeros(Float64,(dim,dim))
  Jac⁻¹ = zeros(Float64,(dim,dim))
  #
  ∇U₁    = zeros(ComplexF64,(dim,dim))
  ∇U₂    = zeros(ComplexF64,(dim,dim))
  sym∇U₁ = zeros(ComplexF64,(dim,dim))
  sym∇U₂ = zeros(ComplexF64,(dim,dim))
  e₁₂    = zeros(ComplexF64,(dim,dim))
  #
  εᵛ₁ = zeros(ComplexF64,dim*(dim-1)) 
  εᵛ₂ = zeros(ComplexF64,dim*(dim-1)) 
  eᵛ₁₂ = zeros(ComplexF64,dim*(dim-1)) 
  #
  σᵛ₁ = zeros(ComplexF64,dim*(dim-1))
  σᵛ₂ = zeros(ComplexF64,dim*(dim-1))
  σᵛ₁₂ = zeros(ComplexF64,dim*(dim-1))
  #
  sym∇ = zeros(ComplexF64,(nne_max*dim,dim*(dim-1)))
  sym∇ⁿˡ₁ = zeros(ComplexF64,(nne_max*dim,dim*(dim-1)))
  sym∇ⁿˡ₂ = zeros(ComplexF64,(nne_max*dim,dim*(dim-1)))
  #
  for iΩ ∈ mesh.Ω
    Dᵢⱼₖₗ = iΩ.mat.Dᵢⱼₖₗ
    for set = 1:iΩ.Sen
      etype = iΩ.Set[set]
      qr = select_quadrature_rule(etype)
      nn = iΩ.Senn[set]
      skip = iΩ.eskip[set]
      for e = 1:iΩ.ne[set]
        conn = @view iΩ.e2n[skip+1+(e-1)*nn:skip+e*nn]
        get_coor!(mesh, conn, X, nn)
        dofs!(U,nn,conn,dofs)
        #
        @inbounds for i = 1:nn*dim
          if (dofs[i]>0)
            Ψ₁ₑ[i] = Ψ₁[dofs[i]]
            Ψ₂ₑ[i] = Ψ₂[dofs[i]]
          else
            Ψ₁ₑ[i] = 0.0
            Ψ₂ₑ[i] = 0.0
          end
        end
        #
        integrate_G!(Fₑ,X,Dᵢⱼₖₗ,Ψ₁ₑ,Ψ₂ₑ,
                     N,∂N∂a,∂N∂x,Jac,Jac⁻¹,
                     ∇U₁,∇U₂,sym∇U₁,sym∇U₂,e₁₂,
                     εᵛ₁,εᵛ₂,eᵛ₁₂,σᵛ₁,σᵛ₂,σᵛ₁₂,
                     sym∇,sym∇ⁿˡ₁,sym∇ⁿˡ₂,
                     nn,etype,qr)
        #
        for i = 1:nn*dim
          if (dofs[i]>0)
            @inbounds Cp.nlr[dofs[i],entry] += Fₑ[i]*mult
          end
        end
        #
      end
    end
  end
  #
  return nothing
  #
end


"""
> assembly_MCK!(mesh::Grid, U::Field, M::SparseMatrixCSC, C::SparseMatrixCSC, K::SparseMatrixCSC, α::Float64, β::Float64)
It assemblies mass, damping, and stiffness matrices
- mesh : Grid 
- U : displacement field
- M : mass matrix
- C : damping matrix
- K : stiffness matrix
- α : Rayleigh damping mass proportional coefficient
- β : Rayleigh damping stiffness proportional coefficient
"""
function assembly_MCK!(mesh::Grid, U::Field, 
                       M::SparseMatrixCSC,C::SparseMatrixCSC,K::SparseMatrixCSC,
                       α::Float64,β::Float64)
  #
  fill!(K.nzval,0.0)
  fill!(M.nzval,0.0)
  #
  X = zeros(Float64,nne_max*dim)
  dofs = zeros(Int64,nne_max*dim)
  Kₑ = zeros(Float64,(nne_max*dim,nne_max*dim))
  Mₑ = zeros(Float64,(nne_max*dim,nne_max*dim))
  #
  N = zeros(Float64,nne_max)
  ∂N∂a = zeros(Float64,(nne_max,dim))
  ∂N∂x = zeros(Float64,(nne_max,dim))
  sym∇ = zeros(Float64,(nne_max*dim,dim*(dim-1)))
  Jac = zeros(Float64,(dim,dim))
  Jac⁻¹ = zeros(Float64,(dim,dim))
  #
  for iΩ ∈ mesh.Ω
    Dᵢⱼₖₗ = iΩ.mat.Dᵢⱼₖₗ
    ρ = iΩ.mat.ρ
    for set = 1:iΩ.Sen
      etype = iΩ.Set[set]
      qr = select_quadrature_rule(etype)
      nn = iΩ.Senn[set]
      skip = iΩ.eskip[set]
      for e = 1:iΩ.ne[set]
        conn = @view iΩ.e2n[skip+1+(e-1)*nn:skip+e*nn]
        get_coor!(mesh, conn, X, nn)
        dofs!(U,nn,conn,dofs)
        integrate_MK!(Mₑ,Kₑ,X,ρ,Dᵢⱼₖₗ,N,∂N∂a,∂N∂x,sym∇,Jac,Jac⁻¹,nn,etype,qr)
        assembly_mat_sym!(K,Kₑ,dofs,nn*dim)
        assembly_mat_sym!(M,Mₑ,dofs,nn*dim)
      end
    end
  end
  #
  for i = 1:size(K.nzval)[1]
    @inbounds C.nzval[i] = α*M.nzval[i] + β*K.nzval[i]
  end
  #
  return nothing
end


"""
> assembly_mat_sym!(K::SparseMatrixCSC,Kₑ::Matrix{Float64},dofs::Array{Int64},rl::Int64)
It assembles elemental matrices into sparse matrices. Is assumes CSC format with lower triangular.
- K : global matrix 
- Kₑ : elemental matrix
- dofs : degrees of freedom 
- rl : row length
"""
function assembly_mat_sym!(K::SparseMatrixCSC,Kₑ::Matrix{Float64},
                           dofs::Array{Int64},rl::Int64)
  for irow = 1:rl
    if (dofs[irow]>0)
      for jcol = 1:rl
        if (dofs[jcol]>0 && dofs[irow]>=dofs[jcol])
          assembly!(K,Kₑ[irow,jcol],dofs[irow],dofs[jcol])
        end
      end
    end
  end
  return nothing
end


"""
> assembly!(K::SparseMatrixCSC, val::Float64, ir::Int64, jc::Int64)
Assembly value in sparse matrix using binary search
- K : global matrix 
- val : value to assemble
- ir : row entry
- jc : column entry
"""
function assembly!(K::SparseMatrixCSC,val::Float64,ir::Int64,jc::Int64)
  @inbounds l = K.colptr[jc]-1
  @inbounds r = K.colptr[jc+1]
  while (true)
    m = (l+r)÷2
    if (K.rowval[m]>ir)
      r = m
    elseif (K.rowval[m]<ir)
      l = m
    else
      @inbounds K.nzval[m] += val
      return nothing
    end
  end
  #
  return nothing
end