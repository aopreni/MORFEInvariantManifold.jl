"""
Overview: global calls to perform model order reduction
"""


"""
> MORFE_mech_autonomous(mesh_file,domains_list,materials,boundaries_list,constrained_dof,bc_vals,α,β,Φₗᵢₛₜ,style,max_order,neig=0,nls=0)
Parametrisation of the autonomous system:

\$\\mathbf{M}\\ddot{\\mathbf{U}}+\\mathbf{C}\\dot{\\mathbf{U}}+\\mathbf{F(\\mathbf{U})}=0\$

- mesh_file : name of the mesh file with extension
- domains_list : tags of domains included in the analysis
- materials : materials associated to each domain
- boundaries_list : list of all boundaries included in the analysis
- constrained_dof : list of constrained degrees of freedom
- bc_vals : value of the boundary condition
- α : mass-proportional damping coefficient
- β : stiffness-proportional damping coefficient
- Φₗᵢₛₜ : list of the master modes included in the analysis
- style : parametrisation style
- max_order : maximum order of the asymptotic development
- neig : number of eigenvalues. If not used, then neig = maximum(Φₗᵢₛₜ)
- nls : nonlinear static analysis. 0 = no, 1 = yes.
"""
function MORFE_mech_autonomous(mesh_file,domains_list,materials,
                               boundaries_list,constrained_dof,bc_vals,
                               α,β,
                               Φₗᵢₛₜ,style,max_order,neig=0,nls=0)
  # initialize output folder
  odir = make_output_dir("tmp")
  odir_C, odir_W, odir_f, odir_M = structure_odir(odir,0)
  # basic data structures
  mesh = read_mesh(mesh_file, domains_list, boundaries_list, materials, constrained_dof, bc_vals)
  n2n = NodalConnectivity(mesh)
  U = Field(mesh,3)
  # matrices
  K = init_symCSC(mesh,U,n2n,0,"r")
  M = deepcopy(K)
  C = deepcopy(K)
  assembly_MCK!(mesh,U,M,C,K,α,β)
  K = Symmetric(K,:L)
  M = Symmetric(M,:L)
  C = Symmetric(C,:L)
  #
  if (neig==0)
    neig = maximum(Φₗᵢₛₜ)
  end
  # compute eigenvalues
  λ, ϕ = eigs(K, M, nev=neig, which=:SM)
  λ = real(λ)
  ϕ = real(ϕ)
  mass_normalization!(ϕ,M,neig)
  # export eigenvalues
  export_eig(mesh, U, λ, ϕ, odir, 0)
  # compute linear system quantities
  nm = size(Φₗᵢₛₜ)[1]
  ndofs = nm*2
  neq = K.data.m
  ω₀ = zeros(Float64,nm)
  ζ₀ = zeros(Float64,nm)
  for i = 1:nm
    ω₀[i] = sqrt(λ[Φₗᵢₛₜ[i]])
    ζ₀[i] = 0.5*(α/ω₀[i]+β*ω₀[i])
  end
  # Preallocare arrays and matrices of interest
  sys_rhs = Array{ComplexF64}(undef,neq+ndofs)
  sys_mat = init_symCSC(mesh,U,n2n,ndofs,"c")
  sys_mat = Symmetric(sys_mat,:L)
  sys_res = Array{ComplexF64}(undef,neq+ndofs)
  # Init parametrisation structrure
  Cp = [Parametrisation() for i in 0:max_order]
  rdyn = init_rdyn(ndofs)
  #
  if (nls>0)
    initialize_parametrisation!(Cp[1],0,ndofs,U.neq)
  end
  # impose identity tangency with linear eigenvalues
  initialize_parametrisation!(Cp[2],1,ndofs,U.neq)
  identity_tangency!(Cp[2],nm,neq,Φₗᵢₛₜ,ϕ,ω₀,ζ₀,U)
  #
  rmat = init_realification_matrix(ndofs)
  Wr = zeros(ComplexF64,(U.neq*2,Cp[2].nc))
  fr = zeros(ComplexF64,(ndofs,Cp[2].nc))
  realification!(Cp[2],Wr,fr,1,ndofs,U.neq,rmat)
  export_data!(Cp[2],Wr,fr,1,ndofs,U,M,ϕ,neig,odir_C,odir_W,odir_f,odir_M)
  #
  p = 1
  fill_rdyn!(rdyn,ndofs,fr,Cp[p+1],p)
  # higher orders
  for p = 2:max_order
    #
    println("order: "*string(p))
    Cic = zeros(Int64,p)
    initialize_parametrisation!(Cp[p+1],p,ndofs,U.neq)
    #
    cc = 1
    recursive_assembly!(Cp,ndofs,p,mesh,U,Cic,cc)
    solve_homological!(Cp,ndofs,p,sys_mat,sys_rhs,sys_res,Cic,style,M,C,K)
    #
    Wr = zeros(ComplexF64,(U.neq*2,Cp[p+1].nc))
    fr = zeros(ComplexF64,(ndofs,Cp[p+1].nc))
    realification!(Cp[p+1],Wr,fr,p,ndofs,U.neq,rmat)
    export_data!(Cp[p+1],Wr,fr,p,ndofs,U,M,ϕ,neig,odir_C,odir_W,odir_f,odir_M)
    #
    fill_rdyn!(rdyn,ndofs,fr,Cp[p+1],p)
    #
  end
  #
  save_matcont_rdyn(rdyn,ndofs,odir)
  #
  return Cp, rdyn
  #
end


"""
> MORFE_mech_autonomous(mesh_file,domains_list,materials,boundaries_list,constrained_dof,bc_vals,Ω_list,κ_modes,κ_list,κ_phase,α,β,Φₗᵢₛₜ,style,max_order_a,max_order_na,neig=0,nls=0)
Parametrisation of the autonomous system:

\$\\mathbf{M}\\ddot{\\mathbf{U}}+\\mathbf{C}\\dot{\\mathbf{U}}+\\mathbf{F(\\mathbf{U})}=\\mathbf{F}(t)\$

- mesh_file : name of the mesh file with extension
- domains_list : tags of domains included in the analysis
- materials : materials associated to each domain
- boundaries_list : list of all boundaries included in the analysis
- constrained_dof : list of constrained degrees of freedom
- bc_vals : value of the boundary condition
- Ω_list : list of forcing frequencies
- κ_modes : load multiplier associated to each forcing term
- κ_list : shape of each forcing
- κ_phase : identification of the forcing phase. "c" = cosine, "s" = sine
- α : mass-proportional damping coefficient
- β : stiffness-proportional damping coefficient
- Φₗᵢₛₜ : list of the master modes included in the analysis
- style : parametrisation style
- max_order : maximum order of the asymptotic development
- neig : number of eigenvalues. If not used, then neig = maximum(Φₗᵢₛₜ)
- nls : nonlinear static analysis. 0 = no, 1 = yes.
"""
function MORFE_mech_nonautonomous(mesh_file,domains_list,materials,
                                  boundaries_list,constrained_dof,bc_vals,
                                  Ω_list,κ_modes,κ_list,κ_phase,
                                  α,β,
                                  Φₗᵢₛₜ,style,max_order_a,max_order_na,
                                  neig=0,nls=0)
  # First we perform a parametrisation of the autonomous system
  #
  # initialize output folder
  println("Initializing directories")
  odir = make_output_dir("tmp")
  odir_C, odir_W, odir_f, odir_M = structure_odir(odir,0)
  # basic data structures
  println("Reading mesh")
  mesh = read_mesh(mesh_file, domains_list, boundaries_list, materials, constrained_dof, bc_vals)
  println("Computing nodal connectivity")
  n2n = NodalConnectivity(mesh)
  U = Field(mesh,3)
  # matrices
  println("initializing M C K")
  K = init_symCSC(mesh,U,n2n,0,"r")
  M = deepcopy(K)
  C = deepcopy(K)
  println("Assemblying M C K")
  assembly_MCK!(mesh,U,M,C,K,α,β)
  K = Symmetric(K,:L)
  M = Symmetric(M,:L)
  C = Symmetric(C,:L)
  #
  if (neig==0)
    neig = maximum(Φₗᵢₛₜ)
  end
  # compute eigenvalues
  λ, ϕ = eigs(K, M, nev=neig, which=:SM)
  λ = real(λ)
  ϕ = real(ϕ)
  mass_normalization!(ϕ,M,neig)
  # export eigenvalues
  export_eig(mesh, U, λ, ϕ, odir, 0)
  # compute linear system quantities
  nm = size(Φₗᵢₛₜ)[1]
  ndofs = nm*2
  neq = K.data.m
  ω₀ = zeros(Float64,nm)
  ζ₀ = zeros(Float64,nm)
  for i = 1:nm
    ω₀[i] = sqrt(λ[Φₗᵢₛₜ[i]])
    ζ₀[i] = 0.5*(α/ω₀[i]+β*ω₀[i])
  end
  # Preallocare arrays and matrices of interest
  sys_rhs = Array{ComplexF64}(undef,neq+ndofs)
  sys_mat = init_symCSC(mesh,U,n2n,ndofs,"c")
  sys_mat = Symmetric(sys_mat,:L)
  sys_res = Array{ComplexF64}(undef,neq+ndofs)
  # Init parametrisation structrure
  Cp = [Parametrisation() for i in 0:max_order_a]
  rdyn = init_rdyn(ndofs)
  if (nls>0)
    initialize_parametrisation!(Cp[1],0,ndofs,U.neq)
  end
  # impose identity tangency with linear eigenvalues
  initialize_parametrisation!(Cp[2],1,ndofs,U.neq)
  identity_tangency!(Cp[2],nm,neq,Φₗᵢₛₜ,ϕ,ω₀,ζ₀,U)
  #
  rmat = init_realification_matrix(ndofs)
  Wr = zeros(ComplexF64,(U.neq*2,Cp[2].nc))
  fr = zeros(ComplexF64,(ndofs,Cp[2].nc))
  realification!(Cp[2],Wr,fr,1,ndofs,U.neq,rmat)
  export_data!(Cp[2],Wr,fr,1,ndofs,U,M,ϕ,neig,odir_C,odir_W,odir_f,odir_M)
  #
  p = 1
  fill_rdyn!(rdyn,ndofs,fr,Cp[p+1],p)
  # higher orders
  for p = 2:max_order_a
    #
    println("Order: "*string(p))
    #
    Cic = zeros(Int64,p)
    initialize_parametrisation!(Cp[p+1],p,ndofs,U.neq,'0')
    #
    cc = 1
    recursive_assembly!(Cp,ndofs,p,mesh,U,Cic,cc)
    solve_homological!(Cp,ndofs,p,sys_mat,sys_rhs,sys_res,Cic,style,M,C,K)
    #
    Wr = zeros(ComplexF64,(U.neq*2,Cp[p+1].nc))
    fr = zeros(ComplexF64,(ndofs,Cp[p+1].nc))
    realification!(Cp[p+1],Wr,fr,p,ndofs,U.neq,rmat)
    export_data!(Cp[p+1],Wr,fr,p,ndofs,U,M,ϕ,neig,odir_C,odir_W,odir_f,odir_M)
    #
    fill_rdyn!(rdyn,ndofs,fr,Cp[p+1],p)
    #
  end
  #
  Wr = nothing
  fr = nothing
  #
  Cp_na = Dict()
  #
  for ith_Ω = 1:size(Ω_list)[1]
    # initialise new output directory
    odir_C_c, odir_W_c, odir_f_c, odir_M_c = structure_odir(odir,ith_Ω,'c')
    odir_C_s, odir_W_s, odir_f_s, odir_M_s = structure_odir(odir,ith_Ω,'s')
    #
    Cp⁺ = [Parametrisation() for i in 0:max_order_na]
    #
    append_rdyn_frequency!(rdyn,ndofs,ith_Ω)
    # 
    Ω = Ω_list[ith_Ω]*ω₀[1]
    η = +im*Ω
    for p = 0:max_order_na
      #
      println("Order: "*string(p))
      #
      initialize_parametrisation!(Cp⁺[p+1],p,ndofs,U.neq,1)
      #
      if (p==0)
        Cic = zeros(Int64,1)
        assembly_forcing_vector!(Cp⁺[1],mesh,U,M,ϕ,κ_modes[ith_Ω],κ_list[ith_Ω],κ_phase[ith_Ω])
      else
        Cic = zeros(Int64,p)
      end
      #
      cc = 1
      recursive_assembly_veps!(Cp,Cp⁺,ndofs,p,mesh,U,Cic,cc)
      solve_homological_veps!(η,Cp,Cp⁺,ndofs,p,sys_mat,sys_rhs,sys_res,style,Cic,M,C,K)
      #
      Wr = zeros(ComplexF64,(U.neq*2,Cp⁺[p+1].nc))
      fr = zeros(ComplexF64,(ndofs,Cp⁺[p+1].nc))
      #
      realification!(Cp⁺[p+1],Wr,fr,p,ndofs,U.neq,rmat)
      #
      export_data_veps!(Cp⁺[p+1],Wr,fr,p,ndofs,U,M,ϕ,neig,
                        odir_C_c, odir_W_c, odir_f_c, odir_M_c,
                        odir_C_s, odir_W_s, odir_f_s, odir_M_s)
      #
      fill_rdyn_veps!(rdyn,ndofs,fr,Cp⁺[p+1],p,ith_Ω)
      #
    end
    #
    Cp_na[ith_Ω] = Cp⁺
    #
  end
  #
  save_matcont_rdyn_nonautonomous(rdyn,ndofs,odir,Ω_list)
  #
  return Cp, Cp_na, rdyn
  #
end