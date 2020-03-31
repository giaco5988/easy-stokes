# Easy-stokes

## Introduction

Easy-stokes have been developed during my PhD Thesis at EPFL. It is a Stokes and Laplace solver for axisymmetric flows based on the Boundary Integral Method (BIM) and Spectral Boundary Integral Method (SBIM). In-depth explanations of the mathematical methods used in Easy-stokes can be found in my PhD Thesis.

```
@techreport{gallino2018droplets,
  title={When droplets deform, break up and propel microswimmers},
  author={Gallino, Giacomo},
  year={2018},
  institution={EPFL}
}
```
Easy-drop offers flexibility in defining various axisymmetric geometries. It has been used to compute the motion of deformable interfaces (droplets and bubbles) and artificial microswimmers.

```
@article{gallino2016stability,
  title={The stability of a rising droplet: an inertialess non-modal growth mechanism},
  author={Gallino, Giacomo and Zhu, Lailai and Gallaire, Fran{\c{c}}ois},
  journal={Journal of Fluid Mechanics},
  volume={786},
  year={2016},
  publisher={Cambridge University Press}
}

@article{gallino2018physics,
  title={Physics of Bubble-Propelled Microrockets},
  author={Gallino, Giacomo and Gallaire, Fran{\c{c}}ois and Lauga, Eric and Michelin, Sebastien},
  journal={Advanced Functional Materials},
  volume={28},
  number={25},
  pages={1800686},
  year={2018},
  publisher={Wiley Online Library}
}
```
The Spectral BIM solver allows the user to carry on stability analysis of deformable interfaces and edge-tracking.

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

## Requirements

Matlab2017 or more recent is recommended, with optimization packages installed.

## Installation

Clone the current repository

```
git clone git@github.com:giaco5988/easy-stokes.git
```

## Tutorials

Go to `easy-drop/tutorials` for an hands-on introduction to the code.

* Simulation of a droplet in a capillary tube (possible to add gravity) `easy-drop/tutorials/tutorial-rising-drop-in-pipe`.