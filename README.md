# MORFEInvariantManifold.jl
This package exploits the direct parametrisation method for invariant manifolds to generate accurate **reduced models** of structures subjected to geometric nonlinearities.

> The packages allows parametrizing also the unforced and undamped version of the model to compute the backbone of the system
One particoular feature of this package is that it embeds its own finite element solver. This implies that upon installation the user can create its personalised database for materials and meshes. Furthermore, we provide output and post-process routines compatible with well established continuation packages as **MatCont** for refined analysis on the resulting reduced models. 

> **Warning**: the package is under development, hence it is subjected to frequent updates

## Do I need this package?

The method implemented in this packages provides accurate reduced models of structures that operate at resonance, hence making it amazing for accelerated computation of periodic orbits of finite element models of structures. Furthermore, numerical coefficients estimated by the presented approach are typically less prone to instabilities compared to techniques serving the same purpose e.g., POD-Galerkin hyper reduction techniques, hence making the computation of stability and bifurcation points more accurate. As a result, the presented method proves amazing for the design of systems that operate at resonance as for instance MEMS micromirrros, gyroscopes, and resonators. If your mechanical structure does not operate at resonance, then this method may serve little use to you.

## Intallation

This package is currently an unregistered Julia package, hence to add it to your packages use the URL of the repository:

`] add https://github.com/aopreni/MORFEInvariantManifold.jl`

We plan to add it to the registered Julia package as soon as possible.

## Citing this work

A reference paper on the developed code is in preparation. For the moment, if you use this package for developing your projects, please cite these works:
* A. Opreni, A. Vizzaccaro, C. Touzé, A. Frangi. [*High order direct parametrisation of invariant manifolds for model order reduction of finite element structures: application to generic forcing terms and parametrically excited systems*](https://www.researchsquare.com/article/rs-1359763/v1), Research Square, 1359763, 2022
* A. Vizzaccaro, A. Opreni, L. Salles, A. Frangi, C. Touzé. [*High order direct parametrisation of invariant manifolds for model order reduction of finite element structures: application to large amplitude vibrations and unconvering of a folding point*](https://arxiv.org/abs/2109.10031). arXiv, 2109.10031, 2021

