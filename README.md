
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


### How it works (over the rationals)

For a system `sys` in variables `vars` and parameters `params`, generically
zero-dimensional in `vars`, the discriminant variety is the union of

- `W_c`: the elimination ideal `<sys, det(Jacobian w.r.t. vars)> ∩ QQ[params]`;
- `W_inf`: for each variable `x` of `vars`, the ideal of the leading coefficients
  (polynomials in `params`) of the elements of a Groebner basis of `<sys>` for the
  block ordering `DRL(vars) > DRL(params)` whose leading monomial in `vars` is a
  pure power of `x`.

Neither needs the reduced Groebner basis of the whole system, which is large
and slow to obtain. Both are computed by the same multi-modular procedure
(`src/modular_block.jl`):

1. Homogenize the system with a new variable `h`, the last one of the parameter
   block: `DRL(vars) > DRL(params, h)`. Setting `h = 1` in a Groebner basis of the
   homogenized system gives a (non-reduced) Groebner basis of `<sys>` for the
   block ordering, so no saturation by `h` is needed.
2. Modulo a prime: Groebner basis of the homogenized system, recording its trace.
   Keep only what is needed, with `h = 1`: the elements free of `vars` (for `W_c`),
   or the leading coefficients of the pure powers (for `W_inf`).
3. Reduced Groebner basis of each kept ideal, a small ideal in the parameters,
   again recording its trace.
4. For further primes, replay both traces, and lift to the rationals only the
   coefficients of the kept part (Chinese remainders and rational reconstruction),
   until one more prime confirms the result. By default the primes are replayed
   four at a time (`batch = 4`, an option of `discriminant_variety`);
   `batch = 1` replays them one after the other.

The generic zero-dimensionality is read on the basis modulo the first prime: every
variable of `vars` has a pure power as a leading monomial, and no element is free
of `vars`. A negative answer is confirmed with a second prime before giving up.
Over finite fields the previous method (Groebner basis over the field) is used.
