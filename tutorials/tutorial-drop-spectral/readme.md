# One droplet using Spectral BIM

## Introduction

This tutorial show how to simulate the motion of a deformable droplet in an ambient fluid due to gravity or an externally imposed flow, as an extensional flow. It shows different methods to solve the governing equations: direct numerical simulation (DNS), Newton Method, Continuation Method, Edge Tracking. Find more info and validations in the reference below.

```
@techreport{gallino2018droplets,
  title={When droplets deform, break up and propel microswimmers},
  author={Gallino, Giacomo},
  year={2018},
  institution={EPFL}
}
```

## Non dimensional parameters

The problem is fully defined with the following non-dimensional parameters.

* Capillary number Ca=U mu_2/gamma, it sets the strenght of the imposed flow
* The viscosity ratio lambda=mu_1/mu_2 between inner and outer fluid.

In case of a droplet rising due to gravity, see the reference below for details regarding governing equations and non-dimensionalization.

```
@article{gallino2016stability,
  title={The stability of a rising droplet: an inertialess non-modal growth mechanism},
  author={Gallino, Giacomo and Zhu, Lailai and Gallaire, Fran{\c{c}}ois},
  journal={Journal of Fluid Mechanics},
  volume={786},
  year={2016},
  publisher={Cambridge University Press}
}
```

In case of a droplet in an extensioanl flow, see the reference below for details regarding governing equations and non-dimensionalization.

```
@article{gallino2018edge,
  title={Edge states control droplet breakup in subcritical extensional flows},
  author={Gallino, Giacomo and Schneider, Tobias M and Gallaire, Fran{\c{c}}ois},
  journal={Physical Review Fluids},
  volume={3},
  number={7},
  pages={073603},
  year={2018},
  publisher={APS}
}
```

## Direct numerical simulation (DNS)

Run the script `tutorial_dns.m` to perform simulations after having set the `REPOSITORY_NAME` variable to `path/to/easy-drop`. By default, it computes the deformation of an initially spherical droplet in an extensional flow when Ca=0.1 and lambda=5. The script is densely commented and describes the main functionalities of the solver.

Run the script `post_processing_dns_spectral.m` to perform post-processing after having run `tutorial_drop_spectral_dns.m ` and set the `REPOSITORY_NAME` variable to `path/to/easy-drop`. The result of `tutorial_drop_spectral_dns.m ` is shown below. It shows the velocity residuals and the converged droplet shape. Low values of the residuals, calculated as velocity normal to the interface in the dropl frame, indicates a converged solution.

![](docs/residuals.png)

## Newton method

The [Newton method](https://en.wikipedia.org/wiki/Newton%27s_method) is probably the most common non time-marching method to solve differential equation. It has two main advantages:

*  It converges to numerical precision fast.
*  It converges to steady solutions, no matter if stable or unstable.

It's main drawback compared to DNS is the need of computing the [Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of the time evolution equation, which might be mathematically involved and computationally expensive.

Run the script `tutorial_newton_method.m` to perform simulations after having set the `REPOSITORY_NAME` variable to `path/to/easy-drop`. By default, it computes the shape in an extensional flow when Ca=0.1 and lambda=5. Below the results, residuals over iterations and the converged shape.

![](docs/newton_plot.png)

## Stability analysis

Once the fixed point has been found with the Newton Method, it is insightful to study the [stability](Stability_theory) of the fixed by computing the eigenvalues of the Jacobian. If positive eigenvalues are found, the fixed point is unstable. Run `` to compute the eigenvalue spectra of the solution from the previous point. Below, see the eigenvalues spectra, since all eigenvalues are negative, the solution is stable.

![](docs/stability.png)

## Continuation method

Based on Newton method, one can perform [continuation method](https://en.wikipedia.org/wiki/Numerical_continuation) in order to compute fixed points of the equations vs changes in parameters. Here we implement [pseudo-arclength continuation](https://en.wikipedia.org/wiki/Numerical_continuation#Pseudo-arclength_continuation) in order to turn atound bifurcation in the bifurcation diagram. Below see the bifurcation diagram for and extensional flow for lambda=5. Black dots indicate stable solutions and red dots unstable ones. They meet in a saddle node bifurcation for Ca~0.1.

![](docs/bifurc.png)

## Edge tracking

Edge tracking is a strategy to track unstable solutions living on the boundaries of the basin of attraction, which separates stable from unstable trajectories of the dynamical system, so called edge states.
