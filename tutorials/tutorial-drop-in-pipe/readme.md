# Tutorial drop in pipe

## Introduction

This tutorial simulate the motion of a deformable droplet moving in a pipe due to an imposed flow and gravity. It is inspired and validated against the paper below.

```
@article{lac2009motion,
  title={Motion of a drop along the centreline of a capillary in a pressure-driven flow},
  author={Lac, Etienne and Sherwood, JD},
  journal={Journal of fluid mechanics},
  volume={640},
  pages={27--54},
  year={2009},
  publisher={Cambridge University Press}
}
```

## Problem definition

![](docs/domain_drop_in_pipe.eps)

Referring to the picture above, the dimensioanl parameters of the problem are:

* The droplet viscosity mu_1.
* The ambient fluid viscosity mu_2.
* The capillary tube radius R.
* The droplet volume 4/3*pi*R^3*aplha^3.
* The average inlet velocity U.
* The interface surface tension gamma.
* The accelleration of gravity g.
* The difference between droplet density and ambient fluid density delta-rho.

## Non dimensional parameters
The problem is fully defined with the following non-dimensional parameters, following the reference scales [L]=R, [V]=gamma/mu_2, [P]=gamma/R.

* Bond number Bo=delta-rho R^2 g/gamma, it sets the strenght of the gravity force.
* Capillary number Ca=U mu_2/gamma, it sets the strenght of the imposed flow
* Droplet size as the ratio \alpha between the radius of the droplet (when spherical) and the radius of the capillary tube.
* The viscosity ratio lambda=mu_1/mu_2 between inner and outer fluid.

## Run simulations

Run the script `tutorial_drop_in_pipe.m` to perform simulations after having set the `REPOSITORY_NAME` variable to `path/to/easy-drop`. By default, the script reproduces the results shown in figure-2 of Lac and Sherwood.

The script is densely commented and describes the main functionalities of the solver.

## Post-processing and results

Run the script `post_processing_drop_pipe.m` to perform post-processing after having run `tutorial_drop_in_pipe.m` and set the `REPOSITORY_NAME` variable to `path/to/easy-drop`.

The result of `tutorial_drop_in_pipe.m`, keeping the default values is shown below. It gives the converged droplet shape as in figure-2 of Lac and Sherwood. Low values of the residuals, calculated as velocity normal to the interface in the dropl frame, indicates a converged solution.

![](docs/residuals_drop_in_pipe.eps)
