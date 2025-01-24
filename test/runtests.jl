using Test
using Nemo
using DiscriminantVariety

@testset "stability_bouzidi_rouillier" begin
    include("../Systems/stability_bouzidi_rouillier.jl")
    
    dv = DiscriminantVariety.discriminant_variety(sys, [x,y], [u1,u2])

    @test dv == [[u1 - u2 - 11], [5*u1 + 3*u2 + 2], [u1 - 5*u2 - 1], [4*u1 - 2*u2 + 3], [7*u1 + u2 - 3], [7*u1 - 3*u2 + 7], [6*u1^2 + 4*u1*u2 + 2*u2^2 - 8*u2 + 1], [6*u1^2 - 6*u1*u2 - 4*u2^2 + 25*u1 + 3*u2 + 11], [1276*u1^6 - 2828*u1^5*u2 - 168*u1^4*u2^2 + 2896*u1^3*u2^3 + 1544*u1^2*u2^4 + 340*u1*u2^5 + 76*u2^6 + 874*u1^5 - 10474*u1^4*u2 - 4984*u1^3*u2^2 - 4300*u1^2*u2^3 - 1866*u1*u2^4 + 14*u2^5 - 72*u1^4 - 6542*u1^3*u2 + 6663*u1^2*u2^2 - 1396*u1*u2^3 - 1053*u2^4 - 239*u1^3 - 2461*u1^2*u2 + 8675*u1*u2^2 + 665*u2^3 + 170*u1^2 - 1834*u1*u2 + 2064*u2^2 + 301*u1 - 557*u2 + 91]]
end

