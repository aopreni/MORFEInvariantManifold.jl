# Overview

This package develops reduced models for continuum structures using the Direct Parametrization mehtod for Invarint Manifolds. The method parametrizes the system motion along an invariant set of the system by introducing a nonlinear change of coordinates between the model coordinates (e.g. nodal displacement and velocity) and the reduced model coordinates. For the autonomous system, the coordinate change is formulated as:

$$\mathbf{U} = \boldsymbol{\Psi}(\mathbf{\mathrm{z}}),$$
$$\mathbf{V} = \boldsymbol{\Upsilon}(\mathbf{\mathrm{z}}),$$

with $\mathbf{U}$ displacement field, $\mathbf{V}$ velocity field, and $\boldsymbol{\Psi}(\mathbf{\mathrm{z}})$ and $\boldsymbol{\Upsilon}(\mathbf{\mathrm{z}})$ nonlinear functions of the *normal coordinates* $\mathbf{\mathrm{z}}$. The method provides the dynamics of the system along the embedding as well. This is often referred in literature as *reduced dynamics*:

$$\dot{\mathbf{\mathrm{z}}} = \mathbf{f}(\mathbf{\mathrm{z}})$$

with $\mathbf{f}(\mathbf{\mathrm{z}})$ polynomial law in $\mathbf{\mathrm{z}}$.

In presence of non-autonomous forcing, then the procedure exploits a perturbaitve approach to compute the whisker attached to the manifold computed from the autonomous problem. Let us consider an harmonic forcing characterised by a finite set of harmonic components $\boldsymbol{\Omega}$. The associated expansion for the non-autonomous problem is:

$$\mathbf{U} = \boldsymbol{\Psi}(\mathbf{\mathrm{z}})+\varepsilon\tilde{\boldsymbol{\Psi}}(\mathbf{\mathrm{z}},\boldsymbol{\Omega},t),$$
$$\mathbf{V} = \boldsymbol{\Upsilon}(\mathbf{\mathrm{z}})+\varepsilon\tilde{\boldsymbol{\Upsilon}}(\mathbf{\mathrm{z}},\boldsymbol{\Omega},t),$$

with $\tilde{\boldsymbol{\Psi}}(\mathbf{\mathrm{z}},\boldsymbol{\Omega},t)$ and $\tilde{\boldsymbol{\Upsilon}}(\mathbf{\mathrm{z}},\boldsymbol{\Omega},t)$ nonautonomous coordinate transformation. Similarly, we can introduce a non-autonomous reduced dynamics as:

$$\dot{\mathbf{\mathrm{z}}} = \mathbf{f}(\mathbf{\mathrm{z}}) + \tilde{\mathbf{f}}(\mathbf{\mathrm{z}},\boldsymbol{\Omega},t),$$

which provides the reduced model we need to solve.

The method is based on the substitution of the nonlinear coordinate change in the original equations of motion to retrieve the so-called *invariance equation*. Its solution allows deriving analytic expressions for mappings and reduced dynamics, the latter corresponding to the reduced model of the problem. The appealing feature of the method is that the dimension of the reduced model is equal to the dimension of the invariant set we parametrise our system motion on, hence it is relatively independent from the dimension of the original finite element model. The result is that this method provides the smaller possible reduced model for resonating structures.