
### Installation

In Julia, run the following

```julia
using Pkg; Pkg.add("DiscriminantVariety")
```

or, in your favorite terminal, get the sources

```
git clone https://github.com/sumiya11/DiscriminantVariety.jl
```

Documentation: https://sumiya11.github.io/DiscriminantVariety.jl/

### Usage example

In Julia, type the following

```julia
using Nemo, DiscriminantVariety

R, (x,y,z,a,b,c) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"])
sys = [a*x^2 + b - 1, y + b*z, y + c*z]
@show discriminant_variety(sys, [x,y,z], [a,b,c])
```


### What is computed

For a system `sys` in variables `vars` and parameters `params`, generically
zero-dimensional in `vars`, the discriminant variety is the union of

- `W_c`: the elimination ideal `<sys, det(Jacobian w.r.t. vars)> ∩ QQ[params]`;
- `W_inf`: for each variable `x` of `vars`, the ideal of the leading coefficients
  (polynomials in `params`) of the elements of a Groebner basis of `<sys>` for the
  block ordering `DRL(vars) > DRL(params)` whose leading monomial in `vars` is a
  pure power of `x`.
