
### Installation

In Julia, run the following

```julia
using Pkg; Pkg.add("DiscriminantVariety")
```

or, in your favorite terminal, get the sources

```
git clone https://github.com/sumiya11/DiscriminantVariety.jl
```

### Usage example

In Julia, type the following

```julia
using Nemo, DiscriminantVariety

R, (x,y,z,a,b,c) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"])
sys = [a*x^2 + b - 1, y + b*z, y + c*z]
@show discriminant_variety(sys, [x,y,z], [a,b,c])
```

