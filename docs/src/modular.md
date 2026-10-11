# Computation over QQ

Over ``\mathbb{Q}``, ``W_c`` and ``W_\infty`` are ideals in the parameters. Neither needs the reduced Gröbner basis of ``I = \langle P_1,\ldots,P_s \rangle``, which is large and expensive to compute: only these ideals are lifted, multi-modularly.

1. The system is homogenized with a new variable ``h``, placed last in the block of the parameters: ``\mathrm{DRL}(X) > \mathrm{DRL}(U, h)``. If ``G`` is a Gröbner basis of the homogenized system, ``G`` with ``h = 1`` is a (non-reduced) Gröbner basis of ``I`` for ``\mathrm{DRL}(X) > \mathrm{DRL}(U)``; no saturation by ``h`` is needed.
2. Modulo a prime ``p``, a Gröbner basis of the homogenized system is computed and its trace recorded (trace 1). Only the part of interest is kept, with ``h = 1``: the elements free of ``X`` for ``W_c``, the leading coefficients of the pure powers of each ``X_i`` for ``W_\infty``.
3. The reduced Gröbner basis of each kept ideal, a small ideal in the parameters, is computed and its trace recorded (trace 2).
4. For further primes, trace 1 and trace 2 are replayed and only the coefficients of the kept ideals are lifted to ``\mathbb{Q}`` (Chinese remaindering and rational reconstruction), until one more prime confirms the result.

The generic zero-dimensionality of the system is read on the basis modulo the first prime; a negative answer is confirmed with a second prime.

The keyword `batch` of [`discriminant_variety`](@ref) is the number of primes for which trace 1 is replayed together, one of `1, 2, 4, 8, 16, 32, 64, 128`; `batch = 1` replays one prime after the other. The result does not depend on it.

```julia
using Nemo, DiscriminantVariety

R, (x, y, a, b) = polynomial_ring(QQ, ["x", "y", "a", "b"])
sys = [QQ(1,7)*a*x^6 + b*y^2 - QQ(1,3), x^2 - QQ(1,5)*a*y - b]
discriminant_variety(sys, [x, y], [a, b]; batch = 1)
```
