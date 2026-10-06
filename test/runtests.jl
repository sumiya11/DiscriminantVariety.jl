using Test
using Nemo
import AbstractAlgebra
using DiscriminantVariety
import Groebner

@testset "Basic" begin
    for interface in [Nemo, AbstractAlgebra]
    for k in [interface.QQ,] #interface.GF(2^30+3)]
    
    R, (x,a) = polynomial_ring(k, ["x","a"])
    dv = DiscriminantVariety.discriminant_variety([a*x^2 + a + 1], [x], [a])
    @test Set(dv) == Set([[a], [a + 1]])
    
    end
    end
end

@testset "stability_bouzidi_rouillier" begin
    include("../Systems/stability_bouzidi_rouillier.jl")
    
    dv = DiscriminantVariety.discriminant_variety(sys, [x,y], [u1,u2])

    @test Set(dv) == Set([[u1 - u2 - 11], [5*u1 + 3*u2 + 2], [u1 - 5*u2 - 1], [4*u1 - 2*u2 + 3], [7*u1 + u2 - 3], [7*u1 - 3*u2 + 7], [6*u1^2 + 4*u1*u2 + 2*u2^2 - 8*u2 + 1], [6*u1^2 - 6*u1*u2 - 4*u2^2 + 25*u1 + 3*u2 + 11], [1276*u1^6 - 2828*u1^5*u2 - 168*u1^4*u2^2 + 2896*u1^3*u2^3 + 1544*u1^2*u2^4 + 340*u1*u2^5 + 76*u2^6 + 874*u1^5 - 10474*u1^4*u2 - 4984*u1^3*u2^2 - 4300*u1^2*u2^3 - 1866*u1*u2^4 + 14*u2^5 - 72*u1^4 - 6542*u1^3*u2 + 6663*u1^2*u2^2 - 1396*u1*u2^3 - 1053*u2^4 - 239*u1^3 - 2461*u1^2*u2 + 8675*u1*u2^2 + 665*u2^3 + 170*u1^2 - 1834*u1*u2 + 2064*u2^2 + 301*u1 - 557*u2 + 91]])
end


@testset "modular elimination" begin
    # The multi-modular elimination (homogenized system, only the elimination
    # ideal lifted) must agree with a Groebner basis over QQ.
    function check(sys, vars)
        J = DiscriminantVariety.jacobian_minors(sys, vars, length(vars))
        full = vcat(sys, J)
        a = DiscriminantVariety.eliminate_modular(full, vars)
        b = DiscriminantVariety.eliminate_groebner(full, vars)
        Set(a) == Set(b)
    end

    R1, (x1,y1,z1,a1,b1,c1) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"], internal_ordering=:degrevlex)
    @test check([a1*x1^2 + b1 - 1, y1 + b1*z1, y1 + c1*z1], [x1,y1,z1])

    R2, (x2,y2,a2,b2) = polynomial_ring(QQ, ["x", "y", "a", "b"], internal_ordering=:degrevlex)
    @test check([a2*x2^6 + b2*y2^2 - 1, x2^2 - a2*y2 - b2], [x2,y2])
    # rational coefficients: needs several primes
    @test check([QQ(1,7)*a2*x2^6 + b2*y2^2 - QQ(1,3), x2^2 - QQ(1,5)*a2*y2 - b2], [x2,y2])

    include("../Systems/stability_bouzidi_rouillier.jl")
    @test check(sys, gens(parent(sys[1]))[1:2])
end


@testset "modular W_inf" begin
    # W_inf lifted multi-modularly must agree with the leading coefficients of
    # the pure powers in a Groebner basis over QQ for DRL(vars) > DRL(params).
    function reference(sys, vars, params)
        ordering = Groebner.DegRevLex(vars) * Groebner.DegRevLex(params)
        gb = Groebner.groebner(sys, ordering=ordering)
        out = []
        for v in vars
            lcs = elem_type(parent(sys[1]))[]
            for f in gb
                lt = Groebner.leading_term(f, ordering=ordering)
                ev = first(exponent_vectors(lt))
                xs = [ev[findfirst(==(w), gens(parent(f)))] for w in vars]
                (xs[findfirst(==(v), vars)] > 0 && count(!iszero, xs) == 1) || continue
                lc = zero(parent(f))
                for (c, e) in zip(coefficients(f), exponent_vectors(f))
                    [e[findfirst(==(w), gens(parent(f)))] for w in vars] == xs || continue
                    lc += c * prod(u^e[findfirst(==(u), gens(parent(f)))] for u in params)
                end
                push!(lcs, lc)
            end
            I = Groebner.groebner(lcs, ordering=Groebner.DegRevLex(params) * Groebner.DegRevLex(vars))
            any(is_constant, I) || push!(out, Set(I))
        end
        Set(out)
    end
    function check(sys, vars, params)
        zerodim, W = DiscriminantVariety.infinity_modular(sys, vars, params)
        zerodim && Set(map(Set, W)) == reference(sys, vars, params)
    end

    R1, (x1,y1,z1,a1,b1,c1) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"], internal_ordering=:degrevlex)
    @test check([a1*x1^2 + b1 - 1, y1 + b1*z1, y1 + c1*z1], [x1,y1,z1], [a1,b1,c1])

    R2, (x2,y2,a2,b2) = polynomial_ring(QQ, ["x", "y", "a", "b"], internal_ordering=:degrevlex)
    @test check([a2*x2^6 + b2*y2^2 - 1, x2^2 - a2*y2 - b2], [x2,y2], [a2,b2])
    @test check([QQ(1,7)*a2*x2^6 + b2*y2^2 - QQ(1,3), x2^2 - QQ(1,5)*a2*y2 - b2], [x2,y2], [a2,b2])

    include("../Systems/stability_bouzidi_rouillier.jl")
    xs = gens(parent(sys[1]))
    @test check(sys, xs[1:2], xs[3:4])
end

@testset "not generically zero-dimensional" begin
    R, (x,y,a,b) = polynomial_ring(QQ, ["x", "y", "a", "b"], internal_ordering=:degrevlex)
    # a relation between the parameters
    @test_throws ErrorException DiscriminantVariety.discriminant_variety([x - a, y - b, a - b], [x,y], [a,b])
    # positive dimension in the variables
    @test_throws ErrorException DiscriminantVariety.discriminant_variety([x*y - a], [x,y], [a,b])
end

@testset "replay one prime at a time or by batches" begin
    R, (x,y,a,b) = polynomial_ring(QQ, ["x", "y", "a", "b"], internal_ordering=:degrevlex)
    # rational coefficients: several primes are needed
    sys = [QQ(1,7)*a*x^6 + b*y^2 - QQ(1,3), x^2 - QQ(1,5)*a*y - b]
    full = vcat(sys, DiscriminantVariety.jacobian_minors(sys, [x,y], 2))
    W_c = [DiscriminantVariety.eliminate_modular(full, [x,y]; batch = k) for k in (1, 2, 4, 8)]
    @test all(==(W_c[1]), W_c)
    W_inf = [DiscriminantVariety.infinity_modular(sys, [x,y], [a,b]; batch = k) for k in (1, 4)]
    @test W_inf[1] == W_inf[2]
    @test DiscriminantVariety.discriminant_variety(sys, [x,y], [a,b]; batch = 1) ==
          DiscriminantVariety.discriminant_variety(sys, [x,y], [a,b])
    @test_throws ArgumentError DiscriminantVariety.eliminate_modular(full, [x,y]; batch = 3)
end
