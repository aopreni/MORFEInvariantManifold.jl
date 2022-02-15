# MORFEInvariantManifold.jl
This package exploits the direct parametrisation method for invariant manifolds to generate accurate **reduced models** of structures subjected to geometric nonlinearities obeying to the following formulation:

$$\mathbf{M}\ddot{\mathbf{U}}+\mathbf{C}\dot{\mathbf{U}}+\mathbf{F(\mathbf{U})}= \mathbf{T}(t),$$
with $\mathbf{M}$ mass matrix, $\mathbf{C}$ damping matrix, $\mathbf{U}$ nodal displacement vector, $\mathbf{F(\mathbf{U})}$ internal force vector, and $\mathbf{T}(t)$ external forcing vector. If a finite elasticity formulation is used, then the internal force vector is a polynomial function of the displacement field and it can be exactly decomposed as:

$$\mathbf{F(\mathbf{U})}=\mathbf{K}\mathbf{U}+\mathbf{G(\mathbf{U},\mathbf{U})}+\mathbf{H(\mathbf{U},\mathbf{U},\mathbf{U})},$$
where $\mathbf{K}$ is the stiffness matrix, $\mathbf{G(\mathbf{U},\mathbf{U})}$ is the quadratic nonlinearity vector, and $\mathbf{H(\mathbf{U},\mathbf{U},\mathbf{U})}$ is the cubic nonlinearity vector. Damping is added to the mode in form of Rayleigh damping:

$$\mathbf{C} = \alpha\mathbf{M} + \beta \mathbf{K},$$
with $\alpha$ and $\beta$ non-negative scalars. 
> The packages allows parametrizing also the unforced and undamped version of the model to compute the backbone of the system
One particoular feature of this package is that it embeds its own finite element solver. This implies that upon installation the user can create its personalised database for materials and meshes. Furthermore, we provide output and post-process routines compatible with well established continuation packages as **MatCont** for refined analysis on the resulting reduced models. 

## Do I need this package?

The method implemented in this packages provides accurate reduced models of structures that operate at resonance, hence making it amazing for accelerated computation of periodic orbits of finite element models of structures. Furthermore, numerical coefficients estimated by the presented approach are typically less prone to instabilities compared to techniques serving the same purpose e.g., POD-Galerkin hyper reduction techniques, hence making the computation of stability and bifurcation points more accurate. As a result, the presented method proves amazing for the design of systems that operate at resonance as for instance MEMS micromirrros, gyroscopes, and resonators. If your mechanical structure does not operate at resonance, then this method may serve little use to you.

## Intallation

This package is currently an unregistered Julia package, hence to add it to your packages simply add it using the URL of the repository:

`] add https://github.com/aopreni/MORFEInvariantManifold.jl`

We plan to add it to the registered Julia package as soon as possible.

## Citing this work

A reference paper on the developed code is in preparation. For the moment, if you use this package for developing your projects, please cite this work:
* A. Vizzaccaro, A. Opreni, L. Salles, A. Frangi, C. Touzé. [*High order direct parametrisation of invariant manifolds for model order reduction of finite element structures: application to large amplitude vibrations and unconvering of a folding point*](https://arxiv.org/abs/2109.10031). arXiv, 2109.10031, 2021

