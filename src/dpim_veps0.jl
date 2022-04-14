"""
Overview: collections of functions required to solve the ε⁰-invariance equation
"""


"""
> solve_homological!(Cp,ndofs,p,sys_mat,sys_rhs,sys_res,ic,style,M,C,K)
It solves al homological equations for a given order of the asymtptotic developmet of the ε⁰-invariance.
- Cp : parametrisation structure
- ndofs : number of degrees of freedom
- p : order of the asymptotic developmet
- sys_mat : matrix used to solve the homological equations
- sys_rhs : right hand side of the homological equations
- sys_res : residual of the homological equations
- ic : index combination
- style : parametrisation style
- M : mass matrix
- C : damping matrix
- K : stiffness matrix
"""
function solve_homological!(Cp,ndofs,p,sys_mat,sys_rhs,sys_res,ic,
                            style,M,C,K)
  #
  resonant_modes = zeros(Bool,ndofs)
  #
  for i = 1:Cp[p+1].nc
    check = Cp[p+1].conj[i]
    if (check>0)
      σ = 0.0+0im
      for j = 1:p
        ic[j] = Cp[p+1].comb[j,i]
        σ += Cp[2].f[ic[j],ic[j]]
      end
      check_resonances!(Cp,σ,ndofs,style,resonant_modes)
      #
      assembly_sys_mat!(Cp,sys_mat.data,M.data,C.data,K.data,σ,resonant_modes,ndofs)
      assembly_sys_rhs!(Cp,sys_rhs,sys_res,M.data,C.data,i,p,ndofs,σ,resonant_modes)
      #
      sys_res = sparse(sys_mat)\sys_rhs
      #
      @inbounds for j = 1:M.data.m
        Cp[p+1].W[j,i] = Cp[p+1].W[j+M.data.m,i] + σ*sys_res[j]
        Cp[p+1].W[j+M.data.m,i] = sys_res[j]
      end
      #
      @inbounds for j = 1:ndofs
        if (resonant_modes[j]>0)
          Cp[p+1].f[j,i] = sys_res[M.data.m+j]
        end
      end
      #
      @inbounds for j = 1:ndofs
        for jj = 1:M.data.m
          Cp[p+1].W[jj,i] += Cp[p+1].f[j,i]*Cp[2].W[jj+M.data.m,j]
        end
      end
      #
    end
  end
  #
  nm = Int(ndofs/2)
  for i = 1:Cp[p+1].nc
    check = Cp[p+1].conj[i]
    tmp = abs(check)
    #
    if (check<0)
      #
      for j = 1:nm
        Cp[p+1].f[j+nm,i] = conj(Cp[p+1].f[j,tmp])
        Cp[p+1].f[j,i] = conj(Cp[p+1].f[j+nm,tmp])
      end
      #
      @inbounds for j = 1:M.data.m*2
        Cp[p+1].W[j,i] = conj(Cp[p+1].W[j,tmp])
      end
      #
    end
    #
  end
  #
  return nothing
  #
end


"""
> assembly_sys_rhs!(Cp,sys_rhs,sys_res,M,C,comb,p,ndofs,σ,resonant_modes)
It assembles the right hand side of a given homological equation of the ε⁰-developmet
- Cp : parametrisation structure
- sys_rhs : right hand side of the homological equations
- sys_res : residual of the homological equations
- M : mass matrix
- C : damping matrix
- comb : reference monomial index
- p : order of the asymptotic development
- ndofs : number of degrees of freedom
- σ : eingevalues λ summation
- resonant_modes : book-keeping of the set R of resonant modes
"""
function assembly_sys_rhs!(Cp,sys_rhs,sys_res,M,C,comb,p,ndofs,σ,resonant_modes)
  #
  fill!(sys_rhs,0.0)
  fill!(sys_res,0.0)
  nm = Int(ndofs/2)
  #
  @inbounds for irow = 1:M.m
    pos1 = M.colptr[irow]
    pos2 = M.colptr[irow+1]-1
    for j = pos1:pos2
      jcol = M.rowval[j]
      sys_rhs[irow] -= M.nzval[j]*Cp[p+1].W[jcol,comb] + (σ*M.nzval[j]+C.nzval[j])*Cp[p+1].W[jcol+M.m,comb]
      sys_res[irow] -= M.nzval[j]*Cp[p+1].W[jcol+M.m,comb]
      if (irow!=jcol)
        sys_rhs[jcol] -= M.nzval[j]*Cp[p+1].W[irow,comb] + (σ*M.nzval[j]+C.nzval[j])*Cp[p+1].W[irow+M.m,comb]
        sys_res[jcol] -= M.nzval[j]*Cp[p+1].W[irow+M.m,comb]
      end
    end
  end
  #
  for i = 1:ndofs
    if (resonant_modes[i]>0)
      if (i<=nm)
        @inbounds for j = 1:M.m
          sys_rhs[M.m+i] += Cp[2].W[M.m+j,i]*sys_res[j]
        end
      else
        @inbounds for j = 1:M.m
          sys_rhs[M.m+i] += Cp[2].W[M.m+j,i-nm]*sys_res[j]
        end
      end
    end
  end
  #
  @inbounds for i = 1:M.m
    sys_rhs[i] -= Cp[p+1].nlr[i,comb]
  end
  #
  return nothing
  #
end


"""
> assembly_sys_mat!(Cp,sys_mat,M,C,K,σ,resonant_modes,ndofs)
it assembles the linear system associated to a homological equation
- Cp : parametrisation structure
- sys_mat : matrix associated to the homological equation
- M : mass matrix
- C : damping matrix
- K : stiffness matrix
- σ : summation of eigenvalues λ
- resonant_modes : book-keeping of the set R of resonant modes
- ndofs : number of degrees of freedom
"""
function assembly_sys_mat!(Cp,sys_mat,M,C,K,σ,resonant_modes,ndofs)
  #
  fill!(sys_mat.nzval,0.0)
  nm = Int(ndofs/2)
  #
  for j = 1:K.m
    pos1 = K.colptr[j]
    pos2 = K.colptr[j+1]-1
    ne = pos2-pos1
    pos1_s = sys_mat.colptr[j]
    @inbounds for i = 0:ne
      sys_mat.nzval[pos1_s+i] = σ^2.0*M.nzval[pos1+i] + σ*C.nzval[pos1+i] + K.nzval[pos1+i]
    end
  end
  #
  for d = 1:ndofs
    if (!resonant_modes[d])
      pos = sys_mat.colptr[K.m+d]
      sys_mat.nzval[pos] = 1.0
    end
  end
  #
  for d = 1:nm
    if (resonant_modes[d])
      if (resonant_modes[d+nm])
        pos = sys_mat.colptr[K.m+d] + 1
        sys_mat.nzval[pos] = 1.0
      end
    end
  end
  #
  for d = 1:ndofs
    if (resonant_modes[d])
      if (d<=nm)
        λ = conj(Cp[2].f[d,d])
      else
        λ = Cp[2].f[d-nm,d-nm]
      end
      vec = (σ-λ)*(Symmetric(M,:L)*Cp[2].W[M.m+1:M.m*2,d])
      #
      @inbounds for j = 1:K.m
        pos = sys_mat.colptr[j+1]-1-ndofs+d
        sys_mat.nzval[pos] = vec[j]
      end
      #
      pos = sys_mat.colptr[K.m+d]
      sys_mat.nzval[pos] = 1.0
      #
    end
  end
  #
  return nothing
  #
end


"""
> check_resonances!(Cp,σ,ndofs,style,resonant_modes)
It checks which mode is resonant. It complies with the chosen style.
- Cp : parametrisation data structure
- σ : summation of eigenvalues λ
- ndofs : number of degrees of freedom
- style : parametrisation style
- resonant_modes : book-keeping of the set R of resonant modes
"""
function check_resonances!(Cp,σ,ndofs,style,resonant_modes)
  #
  fill!(resonant_modes,false)
  #
  if (style=='g')
    fill!(resonant_modes,true)
  else
    for i = 1:ndofs
      λ = Cp[2].f[i,i]
      if (abs(σ-λ)/abs(λ)<=ϵ_tol)
        resonant_modes[i] = true
        if (style=='r')
          if (i <= Int(ndofs/2))
            resonant_modes[i+Int(ndofs/2)] = true
          else
            resonant_modes[i-Int(ndofs/2)] = true
          end
        end
      end
    end
  end
  #
  return nothing
end


"""
> recursive_assembly!(Cp,ndofs,p,mesh,U,ic,cc,nls=0)
It assembles all right hand sides of the homological equations for a given order p of the expansion of the ε⁰ development
- Cp : parametrisation data structure
- ndofs : number of degrees of freedom
- p : order of the asymptotic development
- mesh : mesh data structure
- U : dispalcement field
- ic : index combination of a given monomial
- cc : counter of the asymptotic development
- nls : tag for the nonlinear static analysis. 0 = yes, 1 = yes
"""
function recursive_assembly!(Cp,ndofs,p,mesh,U,ic,cc,nls=0)
  #
  for i = 1:ndofs
    ic[cc] = i
    if (cc<p)
      cc += 1
      recursive_assembly!(Cp,ndofs,p,mesh,U,ic,cc,nls)
      cc -= 1
    else
      #
      pos = find_index_global(ic,ndofs,p)
      
      entry = Cp[p+1].cmap[pos]
      check = Cp[p+1].conj[abs(entry)]
      #
      if (check>0)
        assembly_quadratic_nl!(Cp,entry,ic,ndofs,p,mesh,U,nls)
        assembly_cubic_nl!(Cp,entry,ic,ndofs,p,mesh,U)
        assembly_μ_ν_resid!(Cp,entry,ic,ndofs,p,U.neq)
      end
      #
    end
  end
  #
  return nothing
  #
end


"""
> assembly_μ_ν_resid!(Cp,entry,ic,ndofs,p,neq)
It assembles μ and ν for a given monomial I of order p for the ε⁰ development. 
- Cp : parametrisation data structure
- entry : entry associated to a monomial combination
- ic : monomial index combination
- ndofs : number of degrees of freedom
- p : order of the asymptotic development
- neq : number of equations
"""
function assembly_μ_ν_resid!(Cp,entry,ic,ndofs,p,neq)
  #
  ent = abs(entry)
  #
  for j = 2:p-1
    for k = 0:p-j
      posf = find_index_rdyn(j,k,ndofs,ic)
      cmf = Cp[j+1].cmap[posf]
      if (cmf>0)
        f = @view Cp[j+1].f[:,cmf]
        for s = 1:ndofs
          posw = find_index_map(j,k,s,p,ndofs,ic)
          cmw = Cp[p-j+2].cmap[posw]
          if (cmw>0)
            W = @view Cp[p-j+2].W[:,cmw]
            @inbounds for i = 1:neq*2
              Cp[p+1].W[i,ent] += f[s]*W[i]
            end
          end
        end
      end
    end
  end
  #
  return nothing
  #
end


"""
> assembly_cubic_nl!(Cp,entry,ic,ndofs,p,mesh,U)
It assembles cubic nonlinearities vector for a given monomial I of order p
- Cp : parametrisation data structure
- entry : index entry 
- ic : index combination of the monomial
- ndofs : number or degrees of freedom
- p : order of the asymptotic development
- mesh : grid data structure
- U : displacement field
"""
function assembly_cubic_nl!(Cp,entry,ic,ndofs,p,mesh,U)
  #
  for j = 1:p-1
    for k = 1:p-j-1
      #
      pos1,pos2,pos3 = find_position_h(j,k,p,ndofs,ic)
      #
      cm1 = Cp[j+1].cmap[pos1]
      cm2 = Cp[k+1].cmap[pos2]
      cm3 = Cp[p-j-k+1].cmap[pos3]
      #
      if (cm1>0 && cm2>0 && cm3>0)
        #
        Ψ₁ = @view Cp[j+1].W[U.neq+1:2*U.neq,cm1]
        Ψ₂ = @view Cp[k+1].W[U.neq+1:2*U.neq,cm2]
        Ψ₃ = @view Cp[p-j-k+1].W[U.neq+1:2*U.neq,cm3]
        #
        assembly_H_nl!(Cp[p+1], abs(entry), Ψ₁, Ψ₂, Ψ₃, mesh, U)
        #
      end
      #
    end
  end
  #
  return nothing
  #
end


"""
> assembly_quadratic_nl!(Cp,entry,ic,ndofs,p,mesh,U,nls)
It assembles quadratic nonlinearities vector for a given monomial I of order p
- Cp : parametrisation data structure
- entry : index entry 
- ic : index combination of the monomial
- ndofs : number or degrees of freedom
- p : order of the asymptotic development
- mesh : grid data structure
- U : displacement field
- nls : tag for nonlinear static analysis. 0 = no, 1 = yes
"""
function assembly_quadratic_nl!(Cp,entry,ic,ndofs,p,mesh,U,nls)
  #
  for j = 1:p-1
    pos1, pos2 = find_position_g(j,p,ndofs,ic)
    #
    cm1 = Cp[j+1].cmap[pos1]
    cm2 = Cp[p-j+1].cmap[pos2]
    #
    if (cm1>0 && cm2>0)
      #
      Ψ₁ = @view Cp[j+1].W[U.neq+1:2*U.neq,cm1]
      Ψ₂ = @view Cp[p-j+1].W[U.neq+1:2*U.neq,cm2]
      #
      assembly_G_nl!(Cp[p+1], abs(entry), Ψ₁, Ψ₂, mesh, U)
      if (nls>0)
        assembly_H_nl!(Cp[p+1], abs(entry), Ψ₁, Ψ₂, Cp[1].W, mesh, U, 3.0)
      end
    end
  end
  #
  return nothing
  #
end


