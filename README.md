# PeriodicSchurBifurcationKit

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev) | [![Build status](https://github.com/rveltz/PeriodicSchurBifurcationKit.jl/workflows/CI/badge.svg)](https://github.com/rveltz/PeriodicSchurBifurcationKit.jl/actions) [![codecov](https://codecov.io/gh/bifurcationkit/AsymptoticNumericalMethod.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/bifurcationkit/PeriodicSchurBifurcationKit.jl) |


This package provides methods for computing Floquet coefficients associated to periodic orbits computed with [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl).
The methods are based on a **Periodic Schur Decomposition**. This is state of the art in the field.

This Periodic Schur Decomposition is provided by the Julia package [PeriodicSchurDecompositions.jl](https://github.com/RalphAS/PeriodicSchurDecompositions.jl).

## Support and citation
If you use this package for your work, we ask that you cite the following paper. Open source development as part of academic research strongly depends on this. Please also consider starring this repository if you like our work, this will help us to secure funding in the future. It is referenced on HAL-Inria as follows:

```
@misc{veltz:hal-02902346,
  TITLE = {{BifurcationKit.jl}},
  AUTHOR = {Veltz, Romain},
  URL = {https://hal.archives-ouvertes.fr/hal-02902346},
  INSTITUTION = {{Inria Sophia-Antipolis}},
  YEAR = {2020},
  MONTH = Jul,
  KEYWORDS = {pseudo-arclength-continuation ; periodic-orbits ; floquet ; gpu ; bifurcation-diagram ; deflation ; newton-krylov},
  PDF = {https://hal.archives-ouvertes.fr/hal-02902346/file/354c9fb0d148262405609eed2cb7927818706f1f.tar.gz},
  HAL_ID = {hal-02902346},
  HAL_VERSION = {v1},
}
```

 
 
> See also [PeriodicSchurDecompositions.jl](https://github.com/RalphAS/PeriodicSchurDecompositions.jl) for how to cite it.

## Installation

To install it, please run

`] add https://github.com/bifurcationkit/PeriodicSchurBifurcationKit.jl`
