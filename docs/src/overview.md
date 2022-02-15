# Overview

This package develops reduced models for continuum structures using a technique called Direct Parametrization of Invarint Manifolds. The method parametrizes the system motion along an invariant manifold of the system by introducing a nonlinear change of coordinates between the model coordinates and the reduced model coordinates.

$$\mathbf{U} = \boldsymbol{\Psi}(\mathbf{\mathrm{z}}), \quad \mathbf{V} = \boldsymbol{\Upsilon}(\mathbf{\mathrm{z}}),$$

with $\mathbf{U}$ displacement field, $\mathbf{V}$ velocity field, and $\boldsymbol{\Psi}(\mathbf{\mathrm{z}})$ and $\boldsymbol{\Upsilon}(\mathbf{\mathrm{z}})$ nonlinear function of the *normal coordinates* $\mathbf{\mathrm{z}}$. Furthermore, the method provides the dynamics of the reduced model:

$$\dot{\mathbf{\mathrm{z}}} = \mathbf{f}(\mathbf{\mathrm{z}})$$

with $\mathbf{f}(\mathbf{\mathrm{z}})$ polynomial law in $\mathbf{\mathrm{z}}$. If the original system is non autonomous, then auxiliary variables are added to recast the system as autonomous. 