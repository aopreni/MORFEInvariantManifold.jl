"""
Overview: We here collect all functions and data structures required to manage
solid materials. The package has an internal materials database that is 
updateed by the user to insert his own materials.
"""

"""
> struct material
- ρ : density of the material
- Dᵢⱼₖₗ : fourth order elasticity tensor in Voigt notation
"""
struct material 
  #
  ρ::Float64
  Dᵢⱼₖₗ::Matrix{Float64}
  #
end


"""
> load_material( name::String )
It loads the data structure from the database
- name : name of the material to load
"""
function MORFE_load_material(name::String)
  pth = joinpath(get_materials_path(),name*".jld2")
  data = load(pth)
  return data["material"]
end


"""
> add_material( name::String, ρ::Float64, E::Float64, ν::Float64 )
It adds an isotropic material to the internal database of the package.
- name = name of the material
- ρ = density of the material
- E = Young's modulus 
- ν = Poisson's ratio
**return**  *nothing*
"""
function MORFE_add_material(name::String, ρ::Float64, E::Float64, ν::Float64)
  #
  pth = get_materials_path()
  #
  Dᵢⱼₖₗ = zeros(Float64,(dim*(dim-1),dim*(dim-1)))
  #
  λ = E*ν/((1.0+ν)*(1.0-2.0*ν))
  μ = E/(1+ν)
  #
  Dᵢⱼₖₗ[1,1] = λ+μ
  Dᵢⱼₖₗ[1,2] = λ
  Dᵢⱼₖₗ[2,1] = λ
  Dᵢⱼₖₗ[2,2] = λ+μ
  Dᵢⱼₖₗ[1,3] = λ
  Dᵢⱼₖₗ[3,1] = λ
  Dᵢⱼₖₗ[2,3] = λ
  Dᵢⱼₖₗ[3,2] = λ
  Dᵢⱼₖₗ[3,3] = λ+μ
  Dᵢⱼₖₗ[4,4] = μ/2.0
  Dᵢⱼₖₗ[5,5] = μ/2.0
  Dᵢⱼₖₗ[6,6] = μ/2.0
  #
  mat = material(ρ, Dᵢⱼₖₗ)
  mname = joinpath(pth,name*".jld2")
  save(mname, Dict("material"=>mat))
  #
  return nothing
end



"""
> add_material( name::String, ρ::Float64, Dᵢⱼₖₗ::Matrix{Float64} )
It adds an anisotropic material to the internal database of the package.
- name = name of the material
- ρ = density of the material
- Dᵢⱼₖₗ = fourth order elasticity tensor
**return**  *nothing*
"""
function MORFE_add_material(name::String, ρ::Float64, Dᵢⱼₖₗ::Matrix{Float64})
  #
  pth = get_materials_path()
  #
  mat = material(ρ, Dᵢⱼₖₗ)
  mname = joinpath(pth,name*".jld2")
  save(mname, Dict("material"=>mat))
  #
  return nothing
end


"""
>  get_materials_path()
It returns the path associated to the material database
**return**  pth::String
"""
function get_materials_path()
  tail = length(joinpath("src","MORFEInvariantManifold.jl"))
  pth = chop(pathof(MORFEInvariantManifold), head = 0, tail = tail)*"materials"
  return pth
end

"""
> list_materials()
It lists all materials in the database.
**return**  *nothing*
"""
function MORFE_list_materials()
  pth = get_materials_path()
  tail = length(".jld2")
  files = readdir(pth)
  for f in files
    println(chop(f,head=0,tail=tail))
  end
  return nothing
end


"""
> delete_materials()
It delets a material in the database.
- name = name of the material to delete
**return**  *nothing*
"""
function MORFE_delete_material(name::String)
  pth = joinpath(get_materials_path(),name)*".jld2"
  try
    rm(pth)
  catch
    println("Material not present")
  end
  return nothing
end