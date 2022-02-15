# MORFEInvariantManifold.jl

This package exploits the direct parametrization of invariant manifolds technique to generate accurate **reduced models** of structures subjected to geometric nonlinearities

$$\mathbf{M}\ddot{\mathbf{U}}+\mathbf{C}\dot{\mathbf{U}}+\mathbf{F(\mathbf{U})}= \mathbf{T}(t)$$

with $\mathbf{M}$ mass matrix, $\mathbf{C}$ damping matrix, $\mathbf{U}$ nodal displacement vector, $\mathbf{F(\mathbf{U})}$ internal force vector, and $\mathbf{T}(t)$ external forcing vector. Currently, the method provides Rayleigh damping as dissipation model, hence the matrix $\mathbf{C}$ has the form:

$$\mathbf{C} = \alpha\mathbf{M} + \beta \mathbf{K}$$

with $\alpha$ and $\beta$ non-negative scalars. The internal force term $\mathbf{F(\mathbf{U})}$ is obtained under the assumption of Saint Venant - Kirchhoff constitutive model:

$$\mathbf{F(\mathbf{U})} = \mathbf{K}\mathbf{U} + \mathbf{G}(\mathbf{U},\mathbf{U}) + \mathbf{H}(\mathbf{U},\mathbf{U},\mathbf{U})$$

> The packages allows parametrizing also the unforced and undamped version of the model

One particoular feature of this package is that it embeds its own finite element solver. This implies that upon installation the user can create its personalised database for materials and meshes. Furthermore, we provide a set of functions to exploit well established continuation packages as **BifurcationKit**, **ManLab**, and **MatCont** for refined analysis on the resulting reduced models. 

