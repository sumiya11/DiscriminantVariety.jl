# Discriminant varieties

## Definition

Let ``P_1,\ldots,P_s \in \mathbb{Q}[U_1,\ldots,U_l][X_1,\ldots,X_n]`` define a zero-dimensional system for almost all values of the parameters ``U_1,\ldots,U_l``. The *discriminant variety* of ``P_1,\ldots,P_s``, which can be computed by means of Gröbner bases, is the zero set of polynomial equations in ``U_1,\ldots,U_l`` such that, over any ball of the parameter space that does not meet it, the zero set of the system defines an analytic covering of the ball.

The connected components of its complement are the regions above which the solutions are stable (the implicit function theorem applies on each leaf). In the real case these are the regions where the polynomials defining the discriminant variety have a constant, non-zero sign: the cells of maximal dimension of a *Cylindrical Algebraic Decomposition* (CAD) of the parameter space adapted to these polynomials. The whole parameter space is the union of these cells and of the discriminant variety.

## What is computed

For a system `sys` in the variables `vars` and the parameters `params`, generically zero-dimensional in `vars`, [`discriminant_variety`](@ref) returns the union of

- ``W_c``: the elimination ideal ``\langle P_1,\ldots,P_s, \det J \rangle \cap \mathbb{Q}[U]``, where ``J`` is the Jacobian matrix of the system with respect to ``X_1,\ldots,X_n``;
- ``W_\infty``: for each variable ``X_i``, the ideal of the leading coefficients (polynomials in ``U``) of the elements of a Gröbner basis of ``\langle P_1,\ldots,P_s \rangle`` for the block ordering ``\mathrm{DRL}(X) > \mathrm{DRL}(U)`` whose leading monomial in ``X`` is a pure power of ``X_i``.

The result is a list of components, each a list of polynomials in the parameters. Constant components are removed, and a component given by a single polynomial is split into its irreducible factors. A system that is not generically zero-dimensional in `vars` raises an error.

The coefficients may be in ``\mathbb{Q}`` (Nemo or AbstractAlgebra) or in a prime field (Nemo).

## Example

Example 1 of [this paper by G. Moroz et al.](https://members.loria.fr/GMoroz/assets/LGJMcca09.pdf):

```julia
using Nemo, DiscriminantVariety

R, (x, y, z, a, b, c) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"])
sys = [a*x^2 + b - 1, y + b*z, y + c*z]
discriminant_variety(sys, [x, y, z], [a, b, c])
```

## References

1. Daniel Lazard, Fabrice Rouillier, *Solving parametric polynomial systems*, Journal of Symbolic Computation, 42(6), 2007, 636–667, [doi:10.1016/j.jsc.2007.01.007](https://doi.org/10.1016/j.jsc.2007.01.007).
