"""
Overview: collections of functions required to solve the ε¹-invariance equation
"""

"""
> solve_homological_veps!(η,Cp,Cp⁺,ndofs,p,sys_mat,sys_rhs,sys_res,style,ic,M,C,K)
It solves al homological equations for a given order of the asymtptotic developmet of the ε¹-invariance.
- η : eigenvalue associated to the non-autonomous forcing
- Cp : autonomous parametrisation data structure
- Cp⁺ :  non-auonomous parametrisation data structure
- ndofs : number of degrees of freedom
- p : order of the asymptotic development
- sys_mat : homological equation matrix
- sys_rhs : homological equation right hand side
- sys_res : homological equation residual
- style : parametrisation style
- ic : index combination
- M : mass matrix
- C : damping matrix
- K : stiffness matrix
"""
function solve_homological_veps!(η,Cp,Cp⁺,ndofs,p,sys_mat,sys_rhs,
                                 sys_res,style,ic,M,C,K)
  #
  resonant_modes = zeros(Bool,ndofs)
  nm = Int(ndofs/2)
  #
  for i = 1:Cp⁺[p+1].nc
    #
    σ⁺ = η
    #
    for j = 1:p
      ic[j] = Cp[p+1].comb[j,i]
      σ⁺ += Cp[2].f[ic[j],ic[j]]
    end
    #
    check_resonances!(Cp,σ⁺,ndofs,style,resonant_modes)
    println(resonant_modes)
    assembly_sys_mat!(Cp,sys_mat.data,M.data,C.data,K.data,σ⁺,resonant_modes,ndofs)
    assembly_sys_rhs_veps!(Cp,Cp⁺,sys_rhs,sys_res,M.data,C.data,i,p,ndofs,σ⁺,resonant_modes)
    #
    sys_res = sparse(sys_mat)\sys_rhs
    #
    @inbounds for j = 1:M.data.m
      Cp⁺[p+1].W[j,i] = Cp⁺[p+1].W[j+M.data.m,i] + σ⁺*sys_res[j]
      Cp⁺[p+1].W[j+M.data.m,i] = sys_res[j]
    end
    for j = 1:ndofs
      if (resonant_modes[j]>0)
        Cp⁺[p+1].f[j,i] = sys_res[M.data.m+j]
      end
    end
    #
    for j = 1:ndofs
      @inbounds for jj = 1:M.data.m
        Cp⁺[p+1].W[jj,i] += Cp⁺[p+1].f[j,i]*Cp[2].W[jj+M.data.m,j]
      end
    end
    #
  end
  #
  return nothing
  #
end


"""
> assembly_sys_rhs!(Cp,Cp⁺,sys_rhs,sys_res,M,C,comb,p,ndofs,σ,resonant_modes)
It assembles the right hand side of a given homological equation of the ε⁰-developmet
- Cp : autonomous parametrisation data structure
- Cp⁺ :  non-auonomous parametrisation data structure
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
function assembly_sys_rhs_veps!(Cp,Cp⁺,sys_rhs,sys_res,
                                M,C,comb,p,ndofs,
                                σ,resonant_modes)
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
      #
      sys_rhs[irow] -= M.nzval[j]*Cp⁺[p+1].W[jcol,comb] + 
                       (σ*M.nzval[j]+C.nzval[j])*Cp⁺[p+1].W[jcol+M.m,comb]
      sys_res[irow] -= M.nzval[j]*Cp⁺[p+1].W[jcol+M.m,comb]
      #
      if (irow!=jcol)
        #
        sys_rhs[jcol] -= M.nzval[j]*Cp⁺[p+1].W[irow,comb] + 
                         (σ*M.nzval[j]+C.nzval[j])*Cp⁺[p+1].W[irow+M.m,comb]
        sys_res[jcol] -= M.nzval[j]*Cp⁺[p+1].W[irow+M.m,comb]
        #
      end
    end
  end
  # resonant modes
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
  # elastic nl
  @inbounds for i = 1:M.m
    sys_rhs[i] -= Cp⁺[p+1].nlr[i,comb]
  end
  #
  return nothing
  #
end


"""
> recursive_assembly_veps!(Cp,Cp⁺,ndofs,p,mesh,U,Cic,cc,nls=0)
It assembles right hand sides of all homological equations for a given order p of the expansion.
- Cp : autonomous parametrisation data structure
- Cp⁺ :  non-auonomous parametrisation data structure
- ndofs : number of degrees of freedom
- p : order of the asymptotic development
- mesh : grid data structure
- U : displacement field
- Cic : monomials index combination
- cc : recursion depth
- nls : nonlinear static solution. 0 = no, 1 = yes.
"""
function recursive_assembly_veps!(Cp,Cp⁺,ndofs,p,mesh,U,Cic,cc,nls=0)
  #
  for i = 1:ndofs
    Cic[cc] = i
    if (cc<p)
      cc += 1
      recursive_assembly_veps!(Cp,Cp⁺,ndofs,p,mesh,U,Cic,cc,nls)
      cc -= 1
    else
      #
      pos = find_index_global(Cic,ndofs,p)
      entry = Cp⁺[p+1].cmap[pos]
      #
      assembly_quadratic_nl_veps!(Cp,Cp⁺,entry,Cic,ndofs,p,mesh,U,nls)
      assembly_cubic_nl_veps!(Cp,Cp⁺,entry,Cic,ndofs,p,mesh,U)
      assembly_μ_ν_resid_veps!(Cp,Cp⁺,entry,Cic,ndofs,p,U.neq)
      # 
    end
  end
  #
  return nothing
  #
end


"""
> recursive_assembly_veps!(Cp,Cp⁺,ndofs,p,mesh,U,Cic,cc,nls=0)
It assembles right hand sides of all homological equations for a given order p of the expansion.
- Cp : autonomous parametrisation data structure
- Cp⁺ :  non-auonomous parametrisation data structure
- entry : array entry of a given monomial
- ic : index combination associated to a monomial
- ndofs : number of degrees of freedom
- p : order of the asymptotic development
- neq : number of equations
"""
function assembly_μ_ν_resid_veps!(Cp,Cp⁺,entry,ic,ndofs,p,neq)
  #
  ent = abs(entry)
  #
  for k = 2:p
    for l = 0:p-k
      posf = find_index_rdyn(k,l,ndofs,ic)
      cmf = Cp[k+1].cmap[posf]
      #
      if (cmf>0)
        f = @view Cp[k+1].f[:,cmf]
        #
        for s = 1:ndofs
          posw = find_index_map(k,l,s,p,ndofs,ic)
          cmw = Cp⁺[p-k+2].cmap[posw]
          if (cmw>0)
            W⁺ = Cp⁺[p-k+2].W[:,cmw]
            @inbounds for i = 1:neq*2
              Cp⁺[p+1].W[i,ent] += f[s]*W⁺[i]
            end
          end
        end
      end
    end
  end
  #
  for k = 0:p-1
    for l = 0:p-k
      #
      posf = find_index_rdyn(k,l,ndofs,ic)
      #
      cmf = Cp⁺[k+1].cmap[posf]
      #
      if (cmf>0)
        f⁺ = @view Cp⁺[k+1].f[:,cmf]
        for s = 1:ndofs
          #
          posw = find_index_map(k,l,s,p,ndofs,ic)
          cmw = Cp[p-k+2].cmap[posw]
          #
          if (cmw>0)
            W = @view Cp[p-k+2].W[:,cmw]
            @inbounds for i = 1:neq*2
              Cp⁺[p+1].W[i,ent] += f⁺[s]*W[i]
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
> assembly_cubic_nl_veps!(Cp,Cp⁺,entry,Cic,ndofs,p,mesh,U)
It computes cubic nonlinearity vectors for all monomials of order p of the non-autonomous asymptotic development
- Cp : autonomous parametrisation data structure
- Cp⁺ :  non-auonomous parametrisation data structure
- entry : array entry of a given monomial
- ic : index combination associated to a monomial
- ndofs : number of degrees of freedom
- p : order of the asymptotic development
- mesh : grid data structure
- U : displacement field
"""
function assembly_cubic_nl_veps!(Cp,Cp⁺,entry,ic,ndofs,p,mesh,U)
  #
  for k = 1:p-1
    for l = 1:p-k
      #
      pos1,pos2,pos3 = find_position_h(k,l,p,ndofs,ic)
      #
      cm1 = Cp[k+1].cmap[pos1]
      cm2 = Cp[l+1].cmap[pos2]
      cm3 = Cp⁺[p-k-l+1].cmap[pos3]
      #
      if (cm1>0 && cm2>0 && cm3>0)
        #
        Ψ₁ = @view Cp[k+1].W[U.neq+1:2*U.neq,cm1]
        Ψ₂ = @view Cp[l+1].W[U.neq+1:2*U.neq,cm2]
        Ψ₃ = @view Cp⁺[p-k-l+1].W[U.neq+1:2*U.neq,cm3]
        #
        assembly_H_nl!(Cp⁺[p+1], abs(entry), Ψ₁, Ψ₂, Ψ₃, mesh, U, 3.0)
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
> assembly_cubic_nl_veps!(Cp,Cp⁺,entry,Cic,ndofs,p,mesh,U)
It computes quadratic nonlinearity vectors for all monomials of order p of the non-autonomous asymptotic development
- Cp : autonomous parametrisation data structure
- Cp⁺ :  non-auonomous parametrisation data structure
- entry : array entry of a given monomial
- ic : index combination associated to a monomial
- ndofs : number of degrees of freedom
- p : order of the asymptotic development
- mesh : grid data structure
- U : displacement field
- nls : nonlinear static solution. 0 = no, 1 = yes
"""
function assembly_quadratic_nl_veps!(Cp,Cp⁺,entry,ic,ndofs,p,mesh,U,nls)
  #
  for k = 1:p
    pos1, pos2 = find_position_g(k,p,ndofs,ic)
    #
    cm1 = Cp[k+1].cmap[pos1]
    cm2 = Cp⁺[p-k+1].cmap[pos2]
    #
    #
    if (cm1>0 && cm2>0)
      #
      Ψ₁ = @view Cp[k+1].W[U.neq+1:2*U.neq,cm1]
      Ψ₂ = @view Cp⁺[p-k+1].W[U.neq+1:2*U.neq,cm2]
      #
      assembly_G_nl!(Cp⁺[p+1], abs(entry), Ψ₁, Ψ₂, mesh, U, 2.0)
      #
      if (nls>0)
        assembly_H_nl!(Cp⁺[p+1], abs(entry), Ψ₁, Ψ₂, Cp[1].W, mesh, U, 6.0)
      end
      #
    end
    #
  end
  #
  return nothing
  #
end
