## What is a Discriminant Variety ?

Suppose that $P_1,\ldots,P_s \in \mathbb{Q} [U_1,\ldots,U_l][X_1,\ldots ,X_n]$ for a zero-dimensional system for almost all parameter values (the $U_1,\ldots,U_l$), the authors define the *discriminant variety* of $P_1,\ldots,P_s$, which can be computed by means of Gröbner bases, as the zero set of some polynomial equations in $U_1,\ldots,U_l$ such that over any ball in the parameter's space that does not meet the discriminant variety, the zero set of the system defines an analytic covering of the ball.

Given a discriminant or a discriminant variety, both defined as a union of zeroes of polynomials depending on the parameters, one has to describe the connected components of its complement, say the regions above which the solutions are mathematically stable (the implicit function theorem can be applied on each leaf). In the real case, we are then interested in the regions where the polynomials defining the discriminant variety are of constant and non null sign, which corresponds to the cells of maximum dimension in a *Cylindrical Algebraic Decomposition* (CAD) of the parameters' space adapted to the polynomials defining the discriminant variety.

The entire parameter's space can then be viewed as the union of these cells and
of the discriminant variety.

## Installation

```
using Pkg; Pkg.add("DiscriminantVariety");
```

## References

1. Daniel Lazard, Fabrice Rouillier, Solving parametric polynomial systems, Journal of Symbolic Computation, Volume 42, Issue 6, 2007, Pages 636-667, https://doi.org/10.1016/j.jsc.2007.01.007.

