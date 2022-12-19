"""
Overview: collections of routines and utilities to
perform nonlinear static analyses
"""


"""
> construct_rotation_matrix(vectROT::Vector)
Compute the rotation matrix used for rotation NL analysis
- vectROT : vector of the rotations provided in the launch script
"""
function compute_static_equilibrium!(K, M, C, U, mesh, nls, vector_rotation, out_dir)
  if nls>0
    rotation_matrix = build_rotation_matrix(vector_rotation)
    #
    N = deepcopy(C)
    #
    fROT = zeros(Float64, K.data.m)
    #
    assembly_rotation!(mesh, U, N, fROT, rotation_matrix)
    #
    N = Symmetric(N, :L)
    #
    U0 = - (K + N) \ fROT
    update_field!(U, U0)
    #
    solve_nl_static!(mesh, U, fROT, K, N)
    #
    export_U(mesh, U, out_dir)
  end
  #
  return nothing
  #
end


"""
> add_stat_2_param(Cp[1], U)
Add the static displacement field to parametrisation
"""
function add_stat_2_param(param::Parametrisation, U::Field)
  #
  ind = 1
  #
  for i = 1 : U.nen
    if (U.dof[i] != -1)
      #
      param.W[ind] += U.val[i]
      ind += 1
      #
    end
  end
  #
  return nothing
  #
end


"""
> build_rotation_matrix(vectROT::Vector)
Compute the rotation matrix used for rotation NL analysis
- vectROT : vector of the rotations provided in the launch script
"""
function build_rotation_matrix(vectROT::Vector)
  #
  Rx = zeros(Int, 3, 3)
  #
  Rx[3, 2] = 1.0
  Rx[2, 3] = - 1.0
  #
  Ry = zeros(Int, 3, 3)
  Ry[1, 3] = 1.0
  Ry[3, 1] = - 1.0
  #
  Rz = zeros(Int, 3, 3)
  Rz[1, 2] = - 1.0
  Rz[2, 1] = 1.0
  #
  rotation_matrix = vectROT[1] * Rx + vectROT[2] * Ry + vectROT[3] * Rz
  #
  return rotation_matrix
  #
end


"""
> solve_nl_static!(mesh, U, fSTAT, K, N)
Solve the static analysis to compute static equilibrium
- mesh : mesh data structure
- U : displacement field
- fROT : vector of centrifugal load
- K : stiffness matrix
- N : spin softening matrix
"""
function solve_nl_static!(mesh::Grid,
                          U::Field,
                          fROT::Vector{Float64},
                          K, N)
  #
  ϵ = ϵ_nl
  rhs = zeros(Float64, K.data.m)
  #
  niter = 0
  target_res = 1.0
  #
  println("")
  println("Nonlinear rotation analysis")
  println("Residual :")
  #
  while (true)
    #
    niter += 1
    assembly_SVK!(mesh, U, K, rhs)
    #
    for i = 1 : size(K.data.nzval)[1]
      K.data.nzval[i] += N.data.nzval[i]
    end
    #
    for i = 1 : K.data.m
      rhs[i] -= fROT[i]
    end
    #
    U0 = zeros(Float64, U.neq)
    #
    counter = 1
    #
    for p = 1 : U.nen
      # if (U.dof[p] != -1)
      if (U.dof[p] > 0)
        U0[counter] = U.val[p]
        counter += 1
      end
    end
    #
    rhs -= N * U0
    #
    residual = norm(rhs)
    #
    if (niter == 1)
      target_res = norm(rhs) * ϵ
      if (target_res < ϵ)
        println("The numerical value of the target residual is small.")
        println("The target residual is set to the default value of " * string(ϵ_nl) * ".")
        target_res = ϵ
      end
    end
    #
    if (residual < target_res && niter > 1)
      println(residual)
      break
    end
    println(string(residual))
    dU = K \ rhs
    increment_field!(U, dU)
  end
  #
  println("")
  #
  return nothing
  #
end


"""
> assembly_SVK!(mesh, U, K, rhs)
Assemble the tangent stiffness matrix
- mesh : mesh data structure
- U : displacement field
- K : stiffness matrix
- rhs : residuat vector
"""
function assembly_SVK!(mesh::Grid, U::Field, K, rhs::Vector{Float64})
  #
  fill!(K.data.nzval, 0.0)
  fill!(rhs, 0.0)
  #
  X    = zeros(Float64, nne_max * dim)
  dofs = zeros(Int64, nne_max * dim)
  #
  Kₑ   = zeros(Float64, (nne_max * dim, nne_max * dim))
  Uₑ   = zeros(Float64, (nne_max * dim))
  Fₑ   = zeros(Float64, (nne_max * dim))
  #
  ∂N∂a = zeros(Float64, (nne_max, dim))
  ∂N∂x = zeros(Float64, (nne_max, dim))
  G    = zeros(Float64, (dim, dim))
  Gi   = zeros(Float64, (dim, dim))
  xU   = zeros(Float64, (dim, dim))
  E    = zeros(Float64, dim * (dim - 1))
  S    = zeros(Float64, dim * (dim - 1))
  PSₛ  = zeros(Float64, (dim, dim))
  #
  for iΩ ∈ mesh.Ω
    Dᵢⱼₖₗ = iΩ.mat.Dᵢⱼₖₗ
    # ρ = iΩ.mat.ρ
    for set = 1 : iΩ.Sen
      etype = iΩ.Set[set]
      qr = select_quadrature_rule(etype)
      nn = iΩ.Senn[set]
      #
      NNN  = zeros(Float64, (nn * dim, dim * (dim - 1)))
      NN   = zeros(Float64, (nn * dim, dim, dim))
      #
      skip = iΩ.eskip[set]
      for e = 1 : iΩ.ne[set]
        conn = @view iΩ.e2n[skip + 1 + (e - 1) * nn : skip + e * nn]
        get_coor!(mesh, conn, X, nn)
        dofs!(U, nn, conn, dofs)
        dofs_vals!(U, nn, conn, dofs, Uₑ)
        #
        integrate_SVK!(Kₑ, Fₑ, X, Uₑ, Dᵢⱼₖₗ, ∂N∂a, ∂N∂x, G, Gi, xU, E, S, PSₛ, NNN, NN, etype, nn, qr)
        assembly_mat_sym!(K.data, Kₑ, dofs, nn * dim)
        #
        for i = 1 : nn * dim
          jdof = dofs[i]
          if (jdof>0)
            rhs[jdof] += Fₑ[i]
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
> integrate_SVK!()
Assemble tangent stiffness matrix and associated rhs
"""
function integrate_SVK!(Kₑ::Matrix{Float64}, Fₑ::Vector{Float64}, 
                        X::Array{Float64}, Uₑ::Vector{Float64}, 
                        Dᵢⱼₖₗ::Matrix{Float64}, ∂N∂a::Matrix{Float64}, 
                        ∂N∂x::Matrix{Float64}, G::Matrix{Float64}, 
                        Gi::Matrix{Float64}, xU::Matrix{Float64}, 
                        E::Vector{Float64}, S::Vector{Float64}, 
                        PSₛ::Matrix{Float64}, NNN::Matrix{Float64}, 
                        NN::Array{Float64}, etype::Symbol, nn, qr)
  #
  fill!(Kₑ, 0.0)
  fill!(Fₑ, 0.0)
  #
  for (w, gp) in qr
    #
    ∂N∂a!(∂N∂a,gp,Val{etype})
    #
    fill!(G, 0.0)
    #
    for j = 1 : dim
      for i = 1 : dim
        for k = 1 : nn
          G[i, j] += ∂N∂a[k, j] * X[i + (k - 1) * dim]
        end
      end
    end
    #
    esp  = G[1, 1] * (G[2, 2] * G[3, 3] - G[3, 2] * G[2, 3])
    esp -= G[1, 2] * (G[2, 1] * G[3, 3] - G[3, 1] * G[2, 3])
    esp += G[1, 3] * (G[2, 1] * G[3, 2] - G[3, 1] * G[2, 2])
    #
    Gi[1, 1] = (G[2, 2] * G[3, 3] - G[2, 3] * G[3, 2]) / esp
    Gi[2, 1] = (G[2, 3] * G[3, 1] - G[2, 1] * G[3, 3]) / esp
    Gi[3, 1] = (G[2, 1] * G[3, 2] - G[2, 2] * G[3, 1]) / esp
    Gi[1, 2] = (G[1, 3] * G[3, 2] - G[1, 2] * G[3, 3]) / esp
    Gi[2, 2] = (G[1, 1] * G[3, 3] - G[1, 3] * G[3, 1]) / esp
    Gi[3, 2] = (G[1, 2] * G[3, 1] - G[1, 1] * G[3, 2]) / esp
    Gi[1, 3] = (G[1, 2] * G[2, 3] - G[1, 3] * G[2, 2]) / esp
    Gi[2, 3] = (G[1, 3] * G[2, 1] - G[1, 1] * G[2, 3]) / esp
    Gi[3, 3] = (G[1, 1] * G[2, 2] - G[1, 2] * G[2, 1]) / esp
    #
    fill!(∂N∂x, 0.0)
    #
    for j = 1 : dim
      for k = 1 : dim
        for i = 1 : nn
          ∂N∂x[i, j] += ∂N∂a[i, k] * Gi[k, j]
        end
      end
    end
    #
    fill!(xU, 0.0)
    for i = 1 : dim
      for j = 1 : dim
        for k = 1 : nn
          xU[i, j] += ∂N∂x[k,j] * Uₑ[i + (k - 1) * dim]
        end
      end
    end
    #
    E[1] = xU[1, 1]
    E[2] = xU[2, 2]
    E[3] = xU[3, 3]
    E[4] = xU[1, 2] + xU[2, 1]
    E[5] = xU[2, 3] + xU[3, 2]
    E[6] = xU[1, 3] + xU[3, 1]
    #
    for i = 1 : dim
      E[1] += 0.5 * xU[i, 1] * xU[i, 1]
      E[2] += 0.5 * xU[i, 2] * xU[i, 2]
      E[3] += 0.5 * xU[i, 3] * xU[i,3]
      E[4] += 0.5 * xU[i, 1] * xU[i,2] + 0.5 * xU[i,2] * xU[i,1]
      E[5] += 0.5 * xU[i, 3] * xU[i,2] + 0.5 * xU[i,2] * xU[i,3]
      E[6] += 0.5 * xU[i, 1] * xU[i,3] + 0.5 * xU[i,3] * xU[i,1]
    end
    #
    fill!(NNN, 0.0)
    #
    for i = 1 : nn
      NNN[1 + (i - 1) * dim, 1] += ∂N∂x[i, 1]
      NNN[2 + (i - 1) * dim, 2] += ∂N∂x[i, 2]
      NNN[3 + (i - 1) * dim, 3] += ∂N∂x[i, 3]
      NNN[1 + (i - 1) * dim, 4] += ∂N∂x[i, 2]
      NNN[2 + (i - 1) * dim, 4] += ∂N∂x[i, 1]
      NNN[2 + (i - 1) * dim, 5] += ∂N∂x[i, 3]
      NNN[3 + (i - 1) * dim, 5] += ∂N∂x[i, 2]
      NNN[1 + (i - 1) * dim, 6] += ∂N∂x[i, 3]
      NNN[3 + (i - 1) * dim, 6] += ∂N∂x[i, 1]
      #
      for j = 1 : dim
        NNN[j + (i - 1) * dim, 1] += xU[j, 1] * ∂N∂x[i,1]
        NNN[j + (i - 1) * dim, 2] += xU[j, 2] * ∂N∂x[i,2]
        NNN[j + (i - 1) * dim, 3] += xU[j, 3] * ∂N∂x[i,3]
        NNN[j + (i - 1) * dim, 4] += xU[j, 2] * ∂N∂x[i,1]
        NNN[j + (i - 1) * dim, 4] += xU[j, 1] * ∂N∂x[i,2]
        NNN[j + (i - 1) * dim, 5] += xU[j, 2] * ∂N∂x[i,3]
        NNN[j + (i - 1) * dim, 5] += xU[j, 3] * ∂N∂x[i,2]
        NNN[j + (i - 1) * dim, 6] += xU[j, 1] * ∂N∂x[i,3]
        NNN[j + (i - 1) * dim, 6] += xU[j, 3] * ∂N∂x[i,1]
      end
    end
    #
    fill!(S, 0.0)
    #
    for i = 1 : dim * (dim - 1)
      for j = 1 : dim * (dim - 1)
        S[i] += Dᵢⱼₖₗ[i, j] * E[j]
      end
    end
    #
    PSₛ[1, 1] = S[1]
    PSₛ[2, 2] = S[2]
    PSₛ[3, 3] = S[3]
    PSₛ[1, 2] = S[4]
    PSₛ[2, 1] = S[4]
    PSₛ[2, 3] = S[5]
    PSₛ[3, 2] = S[5]
    PSₛ[1, 3] = S[6]
    PSₛ[3, 1] = S[6]
    #
    for i = 1 : nn
      NN[1 + (i - 1) * dim, 1, 1] = ∂N∂x[i, 1]
      NN[2 + (i - 1) * dim, 2, 2] = ∂N∂x[i, 2]
      NN[3 + (i - 1) * dim, 3, 3] = ∂N∂x[i, 3]
      NN[1 + (i - 1) * dim, 1, 2] = ∂N∂x[i, 2]
      NN[2 + (i - 1) * dim, 2, 3] = ∂N∂x[i, 3]
      NN[3 + (i - 1) * dim, 3, 1] = ∂N∂x[i, 1]
      NN[2 + (i - 1) * dim, 2, 1] = ∂N∂x[i, 1]
      NN[3 + (i - 1) * dim, 3, 2] = ∂N∂x[i, 2]
      NN[1 + (i - 1) * dim, 1, 3] = ∂N∂x[i, 3]
    end
    #
    for i = 1 : nn * dim
      for j = 1 : nn * dim
        #
        for k = 1 : dim
          for l = 1 : dim
            for m = 1 : dim
              Kₑ[i, j] += (NN[i, k, l] * PSₛ[l, m] * NN[j, k, m] ) * w * esp
            end
          end
        end
        #
        for k = 1 : dim * (dim - 1)
          for l = 1 : dim * (dim - 1)
            Kₑ[i, j] += (NNN[i, k] * Dᵢⱼₖₗ[k, l] * NNN[j, l] ) * w * esp
          end
        end
        #
      end
    end
    #
    for i = 1 : nn * dim
      for j = 1 : dim * (dim - 1)
        Fₑ[i] -= (NNN[i, j] * S[j]) * w * esp
      end
    end
    #
  end
  #
  return nothing
  #
end


"""
> export_U(mesh, U, out_dir)
Export static displacement as vtk file
"""
function export_U(mesh::Grid, U::Field, out_dir::String)
  #
  odir = out_dir * "/static_solution"
  #
  try
    mkdir(odir)
  catch
    #
  end
  #
  oeig = odir * "/displacement.vtk"
  #
  io = open(oeig,"w")
  write(io,"# vtk DataFile Version 3.0\n")
  write(io,"eigenfunctions\n")
  write(io,"ASCII\n")
  write(io,"\n")
  write(io,"DATASET UNSTRUCTURED_GRID\n")
  write(io,"POINTS "*string(mesh.nn)*" float\n")
  #
  for n = 1:mesh.nn
    for i = 1:dim
      write(io,string(mesh.n2c[i+(n-1)*dim])*" ")
    end
    write(io,"\n")
  end
  write(io,"\n")
  #
  tcn = 0
  tce = 0
  #
  for i = 1:mesh.nΩ
    for j = 1:mesh.Ω[i].Sen
      nn = mesh.Ω[i].Senn[j]
      tcn += mesh.Ω[i].ne[j]
      tce += mesh.Ω[i].ne[j]*nn
    end
  end
  #
  write(io,"CELLS "*string(tcn)*" "*string(tce+tcn)*"\n")
  #
  for d = 1:mesh.nΩ
    iΩ = mesh.Ω[d]
    for set = 1:iΩ.Sen
      etype = iΩ.Set[set]
      nn = iΩ.Senn[set]
      skip = iΩ.eskip[set]
      #
      print_cell_nodes(io,iΩ,set,nn,skip,Val{etype})
      #
    end
  end
  #
  write(io,"\n")
  write(io,"CELL_TYPES "*string(tcn)*"\n")
  #
  for d = 1:mesh.nΩ
    iΩ = mesh.Ω[d]
    for set = 1:iΩ.Sen
      etype = iΩ.Set[set]
      print_cell_type(io,iΩ,set,Val{etype})

    end
  end
  #
  write(io,"\n")
  write(io,"POINT_DATA "*string(mesh.nn)*"\n")
  #
  write(io,"VECTORS StaticDisplacement float\n")
  #
  for i = 1:mesh.nn
    for j = 1:dim
      write(io,string(U.val[j+(i-1)*dim])*" ")
    end
    write(io,'\n')
  end
  #
  close(io)
  #
  return nothing
  #
end


"""
> assembly_spinsoft_matrix!()
Assemble the spin softening matrix
"""
function assembly_rotation!(mesh::Grid, U::Field, 
                            D::SparseMatrixCSC,f::Vector{Float64},matrix_rotation::Matrix{Float64})
  #
  fill!(D.nzval,0.0)
  fill!(f, 0.0)
  #
  X = zeros(Float64,nne_max*dim)
  dofs = zeros(Int64,nne_max*dim)
  Dₑ = zeros(Float64,(nne_max*dim,nne_max*dim))
  fₑ = zeros(Float64,(nne_max*dim))
  #
  N = zeros(Float64,nne_max)
  ∂N∂a = zeros(Float64,(nne_max,dim))
  ∂N∂x = zeros(Float64,(nne_max,dim))
  Jac = zeros(Float64,(dim,dim))
  Jac⁻¹ = zeros(Float64,(dim,dim))
  #
  for iΩ ∈ mesh.Ω
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
        integrate_rotation!(Dₑ,fₑ, matrix_rotation,X,ρ,N,∂N∂a,∂N∂x,Jac,Jac⁻¹,nn,etype,qr)
        #
        assembly_mat_sym!(D,Dₑ,dofs,nn*dim)
        #
        for i = 1 : nn * dim
          if (dofs[i] > 0)
            @inbounds f[dofs[i]] += fₑ[i]
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
> integrate_spinsoft_matrix()
Integrate the spin softening elementary matrix and centrifugal vector elementary vector
"""
function integrate_rotation!(Dₑ::Matrix{Float64},fₑ::Vector{Float64}, matrix_rotation::Matrix{Float64},
                                    X::Array{Float64},ρ::Float64,
                                    N::Vector{Float64},∂N∂a::Matrix{Float64},∂N∂x::Matrix{Float64},
                                    Jac::Matrix{Float64},
                                    Jac⁻¹::Matrix{Float64},nn::Int64,etype::Symbol,qr)
  #
  @inbounds for j = 1 : nn * dim
              fₑ[j] = 0.0
              for i = 1 : nn * dim
                Dₑ[i, j] = 0.0
              end
            end
  #
  for (w, gp) in qr
    #
    N!(N, gp, Val{etype})
    ∂N∂a!(∂N∂a, gp, Val{etype})
    #
    Jac_det = metric!(X, ∂N∂a, ∂N∂x, nn, Jac, Jac⁻¹)
    #
    γ = ρ * w * Jac_det
    #
    matrot2 = matrix_rotation ^ 2
    for j = 1 : nn
      #
      @inbounds d = N[j] * γ
      for k = 1 : dim
        for l = 1 : dim
          @inbounds e = 0
          for i = 1 : nn
            @inbounds e += N[i] * X[(i - 1) * dim + l]
          end
          @inbounds fₑ[k + (j - 1) * dim] += d * matrot2[k, l] * e
        end
      end
      #
      for i = j : nn
        @inbounds c = N[i] * N[j] * γ
        for k = 1 : dim
          for l = 1 : dim
            @inbounds Dₑ[k + (i - 1) * dim, k + (j - 1) * dim] += c * matrot2[k, l]
          end
        end
      end
      #
    end
    #
  end
  #
  for j = 1 : nn * dim
    for i = j + 1 : nn * dim
      @inbounds Dₑ[j, i] = Dₑ[i, j]
    end
  end
  #
  return nothing
  #
end