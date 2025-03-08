# PeriodicSchurBifurcationKit

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev) | [![Build status](https://github.com/rveltz/PeriodicSchurBifurcationKit.jl/workflows/CI/badge.svg)](https://github.com/rveltz/PeriodicSchurBifurcationKit.jl/actions) [![codecov](https://codecov.io/gh/bifurcationkit/AsymptoticNumericalMethod.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/bifurcationkit/PeriodicSchurBifurcationKit.jl) |


This package provides methods for computing Floquet coefficients associated to periodic orbits computed with [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl).
The methods are based on a **Periodic Schur Decomposition**. This is state of the art in the field.

This Periodic Schur Decomposition is provided by the Julia package [PeriodicSchurDecompositions.jl](https://github.com/RalphAS/PeriodicSchurDecompositions.jl).

## 📚 Support and citation
If you use `BifurcationKit.jl` in your work, we ask that you cite the following paper on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib). Open source development as part of academic research strongly depends on this. Please also consider starring this repository if you like our work, this will help us to secure funding in the future.

 
 
> See also [PeriodicSchurDecompositions.jl](https://github.com/RalphAS/PeriodicSchurDecompositions.jl) for how to cite it.

## 📦 Installation

To install it, please run

`] add https://github.com/bifurcationkit/PeriodicSchurBifurcationKit.jl`
