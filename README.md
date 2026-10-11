# DiscriminantVariety.jl: discriminant varieties of parametric polynomial systems

[![Runtests](https://github.com/sumiya11/DiscriminantVariety.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/sumiya11/DiscriminantVariety.jl/actions/workflows/Runtests.yml)
[![Documentation](https://github.com/sumiya11/DiscriminantVariety.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/sumiya11/DiscriminantVariety.jl/actions/workflows/Documentation.yml)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://sumiya11.github.io/DiscriminantVariety.jl/)

DiscriminantVariety.jl computes a discriminant variety of a parametric polynomial
system that is generically zero-dimensional in its unknowns: a set of polynomial
equations in the parameters outside of which the number of solutions is locally
constant and the solutions depend continuously on the parameters. Over `QQ` the
computation is multi-modular and lifts only the ideals in the parameters; over a
prime field it uses Groebner bases from [Groebner.jl](https://github.com/sumiya11/Groebner.jl).

It is part of the [PACE.jl](https://pace.gitlabpages.inria.fr/pace.jl) ecosystem,
which re-exports `discriminant_variety`.
You can find the documentation [here](https://sumiya11.github.io/DiscriminantVariety.jl/).

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

## Usage

```julia
using Nemo, DiscriminantVariety

R, (x, y, z, a, b, c) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"])
sys = [a*x^2 + b - 1, y + b*z, y + c*z]
@show discriminant_variety(sys, [x, y, z], [a, b, c])
```

The result is a list of components, each a list of polynomials in the parameters.

## What is computed

For a system `sys` in variables `vars` and parameters `params`, generically
zero-dimensional in `vars`, the discriminant variety is the union of

- `W_c`: the elimination ideal `<sys, det(Jacobian w.r.t. vars)> ∩ QQ[params]`;
- `W_inf`: for each variable `x` of `vars`, the ideal of the leading coefficients
  (polynomials in `params`) of the elements of a Groebner basis of `<sys>` for the
  block ordering `DRL(vars) > DRL(params)` whose leading monomial in `vars` is a
  pure power of `x`.

## Tests

```sh
julia --project=. test/runtests.jl
```
