#
function MORFE_backbone_modal(param_dir,out_sol,::Type{Val{:MATCONT}},tag=1)
  #
  p_a = readdlm(param_dir*"/comb_a_0/monomials.txt")
  m_a = readdlm(param_dir*"/manifold_a_0/psi.txt")
  # retrieve general info
  nc = size(p_a)[1]
  ndofs = size(p_a)[2]
  nm = size(m_a)[2]
  # extract periodic orbits
  po_all = matread(out_sol)
  po_all = po_all["x"]
  # retrieve number of periodic orbits
  npo = size(po_all)[2]
  lpo = size(po_all)[1]
  lpo = Int((lpo-2)/ndofs)
  # initialise empty frf
  frf = zeros(Float64,(npo,nm+1))
  # computed amplitude
  amp = zeros(Float64,nm)
  # extracted variables
  vars = zeros(Float64,ndofs)
  # loop over periodic orbits and extract maximum amplitude
  for i = 1:npo
    ω = 2π/po_all[end-1,i]
    frf[i,1] = ω
    for j = 1:lpo
      for k = 1:ndofs
        vars[k] = po_all[k+(j-1)*ndofs,i]
      end
      fill!(amp,0.0)
      for c = 1:nc
        Πᵢ = 1.0
        for d = 1:ndofs
          Πᵢ *= vars[d]^p_a[c,d]
        end
        for m = 1:nm
          amp[m] += m_a[c,m]*Πᵢ
        end
      end
      for m = 1:nm
        if (abs(amp[m])>frf[i,m+1])
          frf[i,m+1] = abs(amp[m])
        end
      end
    end
  end
  writedlm(joinpath(param_dir,"backbone"*string(tag)*".txt"),frf)
  return nothing
end
#
#
#
function MORFE_frc_modal(param_dir,out_sol,nfreq,::Type{Val{:MATCONT}},tag=1)
  # read autonomous mappings
  p_a = readdlm(param_dir*"/comb_a_0/monomials.txt")
  m_a = readdlm(param_dir*"/manifold_a_0/psi.txt")
  # read nonautonomous mappings
  p_na = readdlm(param_dir*"/comb_a_1/monomials.txt")
  m_na_c = Dict()
  m_na_s = Dict()
  for i = 1:nfreq
    m_na_c[i] = readdlm(param_dir*"/manifold_c_"*string(i)*"/psi.txt")
    m_na_s[i] = readdlm(param_dir*"/manifold_s_"*string(i)*"/psi.txt")
  end
  # retrieve general info
  nc = size(p_a)[1]
  nc_na = size(p_na)[1]
  ndofs = size(p_a)[2]
  nm = size(m_a)[2]
  # extract periodic orbits
  po_all = matread(out_sol)
  po_all = po_all["x"]
  # retrieve number of periodic orbits
  npo = size(po_all)[2]
  lpo = size(po_all)[1]
  lpo = Int((lpo-2)/(ndofs+nfreq*2))
  # initialise empty frf
  frf = zeros(Float64,(npo,nm+1))
  # computed amplitude
  amp = zeros(Float64,nm)
  # extracted variables
  vars = zeros(Float64,ndofs)
  # extracted harmonics
  harm = zeros(Float64,(2,nfreq))
  #
  for i = 1:npo
    # contrary to the undamped case the 
    # frequency is a control parameter
    # so we do not need to compute it
    # from the period
    ω = po_all[end,i]
    frf[i,1] = ω
    for j = 1:lpo
      #
      for k = 1:ndofs
        vars[k] = po_all[k+(j-1)*(ndofs+2*nfreq),i]
      end
      #
      for k = 1:2*nfreq
        harm[k] = po_all[ndofs+k+(j-1)*(ndofs+2),i]
      end
      #
      fill!(amp,0.0)
      #
      for c = 1:nc
        #
        for d = 1:ndofs
          Πᵢ *= vars[d]^p_a[c,d]
        end
        #
        for m = 1:nm
          amp[m] += Πᵢ*m_a[c,d]
        end
      end
      #
      for f = 1:nfreq
        for c = 1:nc_na
          Πᵢ = 1.0
          #
          for d = 1:ndofs
            Πᵢ *= vars[d]^p_c[c,d]
          end
          #
          for m = 1:nm
            amp[m] += m_na_c[f][c,m]*Πᵢ
            amp[m] += m_na_s[f][c,m]*Πᵢ
          end
        end
      end
      #
      for m = 1:nm
        if (abs(amp[m])>frf[i,m+1])
          frf[i,m+1] = abs(amp[m])
        end
      end
      #
    end
  end
  #
  writedlm(joinpath(param_dir,"frc"*string(tag)*".txt"),frf)
  #
  return nothing
end
#