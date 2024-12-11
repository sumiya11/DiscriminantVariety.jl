module DiscriminantVariety

using Groebner, Nemo

# When this becomes a bottleneck, it can be implemented in Groebner.jl.
function eliminate(sys, vars)
    @assert !isempty(sys) && !isempty(vars) && allunique(vars)
    @assert all(x -> parent(x) == parent(sys[1]), vars)
    all_vars = gens(parent(sys[1]))
    vars_block_1 = vars
    vars_block_2 = setdiff(all_vars, vars)
    ordering = DegRevLex(vars_block_1) * DegRevLex(vars_block_2)
    gb = groebner(sys, ordering=ordering)
    gen = filter(poly -> all(x -> degree(poly, x) == 0, vars_block_1), gb)
    gen
end

function jacobian(sys, vars)
    @assert !isempty(sys) && allunique(vars)
    @assert all(x -> parent(x) == parent(sys[1]), vars)
    ring = parent(sys[1])
    J = Matrix{elem_type(ring)}(undef, length(sys), length(vars))
    for (i, eq) in enumerate(sys)
        for (j, var) in enumerate(vars)
            J[i, j] = derivative(eq, var)
        end
    end
    J
end

function discriminant_variety(sys, vars, params)
    # Implement me!
end

### Example ###

using Nemo

R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"])
sys = [x^2 + y + z, x*y + z]

@show eliminate(sys, [z])

@show jacobian(sys, [x, z])

end # module DiscriminantVariety
