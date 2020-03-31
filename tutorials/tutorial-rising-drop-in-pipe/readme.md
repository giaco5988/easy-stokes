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

## Non dimensional parameters
The problem is fully defined with the following non-dimensional parameters (We take )

* Bond number Bo, it sets the strenght of the gravity force.
* Capillary number Ca , it sets the strenght of the imposed flow
* Droplet size as the ratio \alpha between the radius of the droplet (when spherical) and the radius of the capillary tube.
* The viscosity ratio between inner and outer fluid.

## Run simulations

Run the script `tutorial_drop_in_pipe.m` to perform simulations after having set the `REPOSITORY_NAME` variable to `path/to/easy-drop`. By default, the script reproduces the results shown in figure-2 of Lac and Sherwood.

The script is densely commented and describes the main functionalities of the solver.

## Post-processing and results

Run the script `post_processing_drop_pipe.m` to perform post-processing after having run `tutorial_drop_in_pipe.m` and set the `REPOSITORY_NAME` variable to `path/to/easy-drop`.

The result of `tutorial_drop_in_pipe.m`, keeping the default values is shown below. It gives the converged droplet shape as in figure-2 of Lac and Sherwood. Low values of the residuals, calculated as velocity normal to th einterface in the dropl frame, indicates a converged solution.

![overview_paper](docs/domain_drop_in_pipe.eps)
![overview_paper](docs/residuals_drop_in_pipe.eps)
