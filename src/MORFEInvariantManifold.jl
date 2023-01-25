module MORFEInvariantManifold

  using JLD2
  using AsterReader
  using SparseArrays
  using FEMQuad
  using Dates
  using FileIO
  using LinearAlgebra
  using Arpack
  using DelimitedFiles
  using MAT
  using MATLAB

  export MORFE_add_material, MORFE_list_materials, MORFE_delete_material, MORFE_load_material
  export MORFE_mech_autonomous, MORFE_mech_nonautonomous
  export MORFE_integrate_rdyn_backbone, MORFE_compute_backbone_modal
  export MORFE_integrate_rdyn_frc, MORFE_compute_frc_modal

  # collection of parametrisation calls
  include("dpim_calls.jl")
  # constants used throughout the code
  include("constants.jl")
  # shape functions of the finite elements routines
  include("shape_functions.jl")
  # materials definition, storage, and management
  include("materials.jl")
  # grid data structure and functionalities
  include("mesh.jl")
  # node-to-node connectivity routines
  include("nodal_connectivity.jl")
  # field data structure
  include("field.jl")
  # sparse matrices creation and management
  include("spmat.jl")
  # parametrisation structure and associated functions
  include("param_struct.jl")
  # assemblage functions
  include("assembler.jl")
  # elemental integration routines
  include("elemental.jl")
  # post process routines
  include("postprocess.jl")
  # look-up functions for efficient indexes combinations access
  include("lookup.jl")
  # realification procedures
  include("realification.jl")
  # parametrisation functions for ε⁰ autonomous development
  include("dpim_veps0.jl")
  # parametrisation frunctions for the ε¹ non-autonomous development
  include("dpim_veps1.jl")
  # integrators post process routines
  include("postprocess_integrator.jl")
  # nonlinear static functions
  include("nl_static.jl")

end # module
