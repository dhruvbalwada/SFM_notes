---
redirect_from:
  - "/files/inertial-ranges-kolmogrov"
title: 'Inertial Ranges'
prev_page:
  url: /intro.html
  title: 'Home'
next_page:
  url: /files/vertical_velocities.html
  title: 'Vertical Velocity'
comment: "***PROGRAMMATICALLY GENERATED, DO NOT EDIT. SEE ORIGINAL FILES IN /content***"
---
## Inertial Range Arguments by Kolmogrov

### What is energy as a function of scale?

*Just some definitions here.*

Energy is

$$ E = \frac{1}{2} \overline{\mathbf{u}.\mathbf{u}} $$,

where the average is an ensemble average, sometimes approximated by a time average.

Using the definition of Fourier transform
$$ \mathbf{u}(\mathbf{x}) = \sum_{\mathbf{k}} \hat{\mathbf{u}}(\mathbf{k}) e^{i \mathbf{k.x}}$$,
we can write
$$ \overline{\mathbf{u}.\mathbf{u}} = \overline{\mathbf{u}.\mathbf{u}^*} =  \overline{\sum_{\mathbf{k}} \hat{\mathbf{u}}(\mathbf{k}) e^{i \mathbf{k.x}}. \sum_{\mathbf{k'}} \hat{\mathbf{u}}^*(\mathbf{k'}) e^{-i \mathbf{k'.x}}}$$

$$ = \sum_{\mathbf{k}} \sum_{\mathbf{k'}} \overline{ \hat{\mathbf{u}}(\mathbf{k}) .  \hat{\mathbf{u}}^*(\mathbf{k'})} e^{i (\mathbf{k} - \mathbf{k'}).\mathbf{x}} $$

Assuming *homogeneity*, the LHS $\overline{\mathbf{u}.\mathbf{u}}$ can not be a function of $\mathbf{x}$. This implies that the RHS can not be a function of $\mathbf{x}$ either, and

$$ \overline{\hat{\mathbf{u}}(\mathbf{k}) .  \hat{\mathbf{u}}^* (\mathbf{k'})} = \delta_{\mathbf{kk'}} |\hat{\mathbf{u}}(\mathbf{k})|^2  = 2 \delta_{\mathbf{kk'}} E( \mathbf{k}) \geq 0  $$

$$ E = \frac{1}{2} \overline{\mathbf{u}.\mathbf{u}}  = \sum_{\mathbf{k}} E( \mathbf{k}) $$

Note: This is similar to Parseval's theorem, but the assumption of homogeneity allowed us to get away without doing a spatial average.

Now, assuming *isotropy* implies that $E( \mathbf{k})$ only depends on the magnitude of $\kappa = \sqrt{k^2 + l^2 + m^2}$, and doesn't care about the angle of the wavenumber vector. So,
$$E = \int_0^{\infty} E_o(\kappa) d\kappa = \int \int \int E(\mathbf{k}) dk dl dm )  = 4 \pi \int_0^{\infty} E'(\kappa) \kappa^2 d\kappa$$,
where
$$ E_o(\kappa) = 4 \pi \kappa^2 E'(\kappa)$$ .

Now that we have the definitions out of the way. Let's proceed.

### Kolmogrov's Inertial Range 'Theory'

Kolmogrov considered turbulence in a very large domain, such that the influence of the boundaries was negligible. There is a length scale, $L_F$, at which energy being input - for example a stirrer stirring fluid at a single wavelength and at a rate $\epsilon$. Also, there is a second, and smaller, length scale, $l_\nu$, where energy is being lost to dissipation.

For this simple setup parametric insight says that $E_o(\kappa)$ will be a function of 5 parameters,
$$ E_o (\kappa) = f(\epsilon, \kappa, L_F, l_\nu, \nu) .$$
Further, if we limit ourselves to a range of length scales that is much smaller than the forcing scale and much larger than the damping scale,
$$ E_o (\kappa) = f_I(\epsilon, \kappa) .$$

Now $[\epsilon] = L^2/T^3$, $[\kappa]=1/L$, and $[E_o(\kappa)] = L^3/T^2$, which gives
$$ E_o (\kappa) = C \epsilon^{2/3}\kappa^{-5/3},$$
Kolmogrov's celebrated 5/3 law.

Using the same analysis as above, we can also guesstimate the length scale at which dissipation becomes important
$$ l_\nu = g(\nu, \epsilon) = \frac{\nu^{3/4}}{\epsilon^{1/4}}.$$

Note that the higher the energy injection rate, the smaller is this dissipation length scale, but this variation is very slowly.

It is instructive to look at the total energy in the inertial range,
$$\int_{2\pi/L_F}^{2\pi/l_\nu} E_o(\kappa) d\kappa  = C \epsilon^{2/3} \int_{2\pi/L_F}^{2\pi/l_\nu} \kappa^{-5/3} d \kappa $$
$$ = \frac{3}{2}C \epsilon^{2/3} (L_F^{2/3} - l_\nu^{2/3}).$$
Thus, the total energy containted in the inertial range is finite even as $nu \rightarrow 0$.

However, this is not the case for the enstrophy spectrum $ \vert
{\mathbf{\omega}}\vert ^2 = \kappa^2 = C \epsilon^{2/3}\kappa^{1/3} $. This, function does not have a finite limit as viscosity goes to zero, and the total enstrophy increases unboundedly as Reyonld's number increases.
