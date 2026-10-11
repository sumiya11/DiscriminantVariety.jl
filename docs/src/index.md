# DiscriminantVariety.jl: discriminant varieties of parametric polynomial systems

`DiscriminantVariety.jl` computes a discriminant variety of a parametric polynomial system that is generically zero-dimensional in its unknowns.

It provides:

- the [discriminant variety](discriminant.md) of a system in ``\mathbb{Q}[U_1,\ldots,U_l][X_1,\ldots,X_n]`` or over a prime field, with its two parts: critical values and values at infinity;
- a [multi-modular computation over ``\mathbb{Q}``](modular.md) that lifts only the ideals in the parameters, not the Groebner basis of the system.

It is part of the [PACE.jl](https://pace.gitlabpages.inria.fr/pace.jl) ecosystem, which re-exports [`discriminant_variety`](@ref).

## Installation

DiscriminantVariety is in the General registry:

```julia
using Pkg
Pkg.add("DiscriminantVariety")
```

Then load it:

```julia
using DiscriminantVariety
```
