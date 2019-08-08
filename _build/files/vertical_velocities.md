---
redirect_from:
  - "/files/vertical-velocities"
title: 'Vertical Velocity'
prev_page:
  url: /files/inertial_ranges_Kolmogrov.html
  title: 'Inertial Ranges'
next_page:
  url: 
  title: ''
comment: "***PROGRAMMATICALLY GENERATED, DO NOT EDIT. SEE ORIGINAL FILES IN /content***"
---
## Vertical Velocities

The vertical component of the velocity has received a lot of attention$^1$ in the oceanographic literature in the last few years. Here we show a systematic way to decompose the vertical velocity into parts arising from the velocity that slides along isopycnals, and from the velocity that pushes across isopycnals.

To do this, we need to decompose the full velocity vector into the along and across isopycnal parts. Here we follow Young (2012) to do this projection, as it is quite systematic and simple.

Firstly we define the two coordinate systems; the cartesian coordinate $(x,y,z,t)$ with unit vectors $(\mathbf{i},\mathbf{j},\mathbf{k})$, and the buoyancy coordinate $(\tilde{x}, \tilde{y}, \tilde{b},\tilde{t})$. The non-orthogonal buoyancy coordinate system has two sets of basis vectors: $(\mathbf{e}^1,\mathbf{e}^2,\mathbf{e}^3)=(\mathbf{i},\mathbf{j},\mathbf{\nabla}b)$, and its dual $(\mathbf{e}_1,\mathbf{e}_2,\mathbf{e}_3)= (\mathbf{i}+ \zeta_{\tilde{x}}\mathbf{k} ,\mathbf{j}+ \zeta_{\tilde{y}}\mathbf{k}, \sigma \mathbf{k})$, where $\zeta_{\tilde{i}}$ is the $i$ derivative of the functional representation of the vertical height of an isopycnal and, $\sigma = \zeta_{\tilde{b}}=1/b_z$, is the isopycnal thickness (or the Jacobian of the transformation). Note that $\mathbf{e}_1, \mathbf{e}_2$ are tangent to the isopycnal, while $\mathbf{e}^3$ is normal to the buoyancy surface.

Any vector $\mathbf{u}$, which is independent of the coordinate, can be represented in any of these three basis vector sets, and converted using the bi-orthogonality condition $(\mathbf{e}^i.\mathbf{e}_j = \delta_{ij})$, where $\delta_{ij}$ is the Kronecker delta.

In cartesian coordinates, the velocity vector can be written as
$$ \mathbf{u} = u\mathbf{i}+v\mathbf{j}+w\mathbf{k}$$,
where the vertical velocity can be written using the buoyancy equation as
$$ w = -1/b_z( b_t + ub_x + vb_y  - \overline{\omega} ) = \zeta_{\tilde{t}} + u \zeta_{\tilde{x}} + v \zeta_{\tilde{y}} + \overline{\omega}\zeta_{\tilde{b}}$$
(follow derivation of equation 33 in Young 2012). The representation of the vertical velocity in terms of $\zeta$ makes it clear that the vertical velocity results from local transience, horizontal advection and diabatic diffusion of an isopycnal.

To get the components of velocity in the directions tangent to the isopycnal, we calculate the contravariant components
$$\mathbf{u} = u \mathbf{e}_1 + v \mathbf{e}_2 + \sigma^{-1}(\zeta_{\tilde{t}} + \overline{\omega}\zeta_{\tilde{b}}) \mathbf{e}_3. $$
This decomposition separates the part of the velocity that is along the isopycnal ("lies on the isopycnal"), and the part that is a result of transience and diabatic effects. $u$ and $v$ are the magnitudes of the horizontal components of the velocity in cartesian coordinates, and also the magnitudes of the along isopycnal components of the velocity in buoyancy coordinate. This confusing fact can be understood as the $\mathbf{e}_1$ and $\mathbf{e}_2$ are not of unit magnitude. Also note that $\mathbf{e}_3$ points vertically up, and not perpendicular the buoyancy surface.

So,
$$\mathbf{u}_{along-iso} = u \mathbf{i} + v \mathbf{j} + (u,v).(\zeta_{\tilde{x}},\zeta_{\tilde{y}}) \mathbf{k}  ,$$
which directly indicated the part of the vertical velocity that is a result of sliding along isopycnals. And the remaining part of the velocity, the $\mathbf{e}_3$ component, has to do with diffusion and transience.

Some differences (To-Do):
- This treatement is different from that present in Freilich and Mahadevan 2019. Not convinced if their treatement is correct.
- Confused about how waves move stuff, vs move shapes. Is all of wave phase propogration captures in $\mathbf{e}_3$, or is part of it in what I call $\mathbf{u}_{along-iso}$ too? Had started working on an example with internal waves (following Wolfe 2014) --- important to distinguish between phase motion and particle motion, part of the particle motion is along isopycnals and part of it moves the isopycnals.
- Can we connect the matrix based rotation (Shafer's notes) to the vector based rotations? (Proper treatement is in Griffies 2014, but will need some time investement.)



*Footnote: $^1$maybe unnecessary (we should be thinking of the velocity as a whole, specially given the ocean's strong adiabatic constrain),*
