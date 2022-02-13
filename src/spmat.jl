"""
Overview: functions used to manage sparse matrices
in the code.
"""

"""
> init_symCSC(mesh::Grid, ϕ::Field, n2n::NodalConnectivity, naux::Int64, tp::String)
It initializes a symmetric CSC matrix from a field. It eventually add auxiliary equations.
- mesh : Grid
- ϕ : field of interest 
- nn : number of queries
- n2n : node-to-node connectivity
- naux : number of auxiliary equations
"""
function init_symCSC(mesh::Grid, ϕ::Field, n2n::NodalConnectivity,
                     naux::Int64, tp::String)
  #
  neq = ϕ.neq + naux
  #
  dofs1 = zeros(Int64,ϕ.cpn)
  dofs2 = zeros(Int64,ϕ.cpn)
  ent = zeros(Int64,neq)
  # count nonzero entries
  nnz = 0
  for n = 1:mesh.nn
    dofs!(ϕ,1,[n],dofs1)
    conn = get_node_connectivity(n2n,n)
    lreg = size(conn)[1]
    for i = 1:ϕ.cpn
      if (dofs1[i]>0)
        for j = 1:ϕ.cpn
          if (dofs1[j]>0 && dofs1[i]>=dofs1[j])
            ent[dofs1[j]] += 1
            nnz += 1
          end
        end
        for j = 1:lreg
          dofs!(ϕ,1,[conn[j]],dofs2)
          for k = 1:ϕ.cpn
            if (dofs2[k]>0 && dofs1[i]>=dofs2[k])
              ent[dofs2[k]] += 1
              nnz += 1
            end
          end
        end
      end
    end
  end
  # add auxiliary
  for i = 1:ϕ.neq
    ent[i] += naux
  end
  for i = 1:naux
    nnz += ϕ.neq+i
    ent[ϕ.neq+i] = naux+1-i
  end
  # allocate arrays
  jcol = zeros(Int64,neq+1)
  irow = zeros(Int64,nnz)
  if (tp=="r")
    val = zeros(Float64,nnz)
  elseif (tp=="c")
    val = zeros(ComplexF64,nnz)
  end
  # fill matrices
  j = 1
  jcol[j] = 1
  for i = 1:neq
    j += ent[i]
    jcol[i+1] = j
  end
  #
  fill!(ent,0.0)
  for n = 1:mesh.nn
    dofs!(ϕ,1,[n],dofs1)
    conn = get_node_connectivity(n2n,n)
    lreg = size(conn)[1]
    for i = 1:ϕ.cpn
      idof = dofs1[i]
      if (idof>0)
        for j = 1:ϕ.cpn
          jdof = dofs1[j]
          if (jdof>0 && dofs1[i]>=dofs1[j])
            iv = ent[jdof] + jcol[jdof]
            ent[jdof] += 1
            irow[iv] = idof
          end
        end
        for j = 1:lreg
          dofs!(ϕ,1,[conn[j]],dofs2)
          for k = 1:ϕ.cpn
            jdof = dofs2[k]
            if (jdof>0 && dofs1[i]>=dofs2[k])
              iv = ent[jdof] + jcol[jdof]
              ent[jdof] += 1
              irow[iv] = idof
            end
          end
        end
      end
    end
  end
  # auxiliary
  for i = 1:ϕ.neq
    iv = ent[i] + jcol[i] - 1
    ent[i] += naux
    for j = 1:naux
      irow[iv+j] = ϕ.neq+j
    end
  end
  for i = 1:naux
    iv = ent[ϕ.neq+i] + jcol[ϕ.neq+i] - 1
    ent[i] += naux+1-i
    for j = i:naux
      irow[iv+j-i+1] = ϕ.neq+j
    end
  end
  #
  for i = 1:neq
    iv1 = jcol[i]
    iv2 = jcol[i+1]-1
    sort!(irow[iv1:iv2])
  end
  #
  return SparseMatrixCSC(neq,neq,jcol,irow,val)
  #
end