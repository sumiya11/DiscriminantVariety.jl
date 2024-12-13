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
    J = Nemo.matrix(ring, J)
    J
end

function jacobian_minors(sys, vars, order)
    J = jacobian(sys, vars)
    Nemo.minors(J, order)
end

function is_univariate_term(term, main_var, vars)
    other_vars = filter(x -> x != main_var, vars)
    degree(term, main_var) > 0 && all(x -> degree(term, x) == 0, other_vars)
end

function remove_vars(term, vars)
    for var in vars
        deg = degree(term, var)
        term = term / var^deg
    end
    term
end

# A hack to reconcile Groebner and Nemo orderings
function first_parametric_term_drl(f, vars, params)
    @assert internal_ordering(parent(f)) == :degrevlex
    ts = collect(terms(f))
    ts = sort(ts, by=t -> monomial(remove_vars(t, params), 1), rev=true)
    ts[1]
end

function is_generically_zerodim(sys, vars, params)
    @assert !isempty(sys) && !isempty(vars) && allunique(vars)
    @assert all(x -> parent(x) == parent(sys[1]), vars)
    @assert all(x -> parent(x) == parent(sys[1]), params)
    @assert internal_ordering(parent(sys[1])) == :degrevlex
    ordering = DegRevLex(vars) * DegRevLex(params)
    gb = groebner(sys, ordering=ordering)
    staircase = map(f -> first_parametric_term_drl(f, vars, params), gb)
    closed_staircase = all(var -> any(lead -> is_univariate_term(lead, var, vars), staircase), vars)
    params_relations = filter(poly -> all(x -> degree(poly, x) == 0, vars), gb)
    zerodim = closed_staircase && isempty(params_relations)
    zerodim
end

function discriminant_variety_generically_zerodim(sys, vars, params)
    @assert internal_ordering(parent(sys[1])) == :degrevlex
    @assert is_generically_zerodim(sys, vars, params)
    
    # The only minor is the determinant
    minors = jacobian_minors(sys, vars, length(vars))
    W_c = [eliminate(vcat(sys, minors), vars)]
    
    ordering = DegRevLex(vars) * DegRevLex(params)
    gb = groebner(sys, ordering=ordering)
    staircase = map(f -> first_parametric_term_drl(f, vars, params), gb)
    lead_univariate = filter(t -> any(var -> is_univariate_term(t, var, vars), vars), staircase)
    lead_coeffs = map(t -> remove_vars(t, vars), lead_univariate)
    W_infty = map(l -> [l], lead_coeffs)

    W_d = vcat(W_c, W_infty)
    
    W_d
end

function postprocess(W_d)
    # Make it look more nice.
    W_d = filter(component -> !any(f -> is_constant(f), component), W_d)
    W_d_ = empty(W_d)
    for component in W_d
        if length(component) == 1
            factors = collect(factor(component[1]))
            for (factor, deg) in factors
                push!(W_d_, [factor^deg])
            end
        else
            push!(W_d_, component)
        end
    end
    W_d = W_d_
    W_d
end

function discriminant_variety(sys, vars, params)
    # Convert to drl
    ring = parent(sys[1])
    K = base_ring(ring)
    ring_drl, xs = polynomial_ring(K, symbols(ring), internal_ordering=:degrevlex)
    sys = map(f -> change_base_ring(K, f, parent=ring_drl), sys)
    vars = map(f -> change_base_ring(K, f, parent=ring_drl), vars)
    params = map(f -> change_base_ring(K, f, parent=ring_drl), params)

    if is_generically_zerodim(sys, vars, params)
        W_d = discriminant_variety_generically_zerodim(sys, vars, params)
    else
        error("Non-zerodim systems are not supported")
    end

    W_d = postprocess(W_d)
    
    W_d
end

export discriminant_variety

### Example ###

# Example 1 and Example 2 from
# https://members.loria.fr/GMoroz/assets/LGJMcca09.pdf

using Nemo

R, (x,y,z,a,b,c) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"])
sys = [a*x^2 + b - 1, y + b*z, y + c*z]
@show discriminant_variety(sys, [x,y,z], [a,b,c])

R, (x,y,a,b) = polynomial_ring(QQ, ["x", "y", "a", "b"])
sys = [a*x^6 + b*y^2 - 1, x^2 - a*y - b]
@show discriminant_variety(sys, [x,y], [a,b])

end # module DiscriminantVariety
