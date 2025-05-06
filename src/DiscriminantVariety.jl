module DiscriminantVariety

import AbstractAlgebra
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
    @assert length(term) == 1
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

function is_generically_zerodim(sys, vars, params)
    @assert !isempty(sys) && !isempty(vars) && allunique(vars)
    @assert all(x -> parent(x) == parent(sys[1]), vars)
    @assert all(x -> parent(x) == parent(sys[1]), params)
    @assert internal_ordering(parent(sys[1])) == :degrevlex
    ordering = DegRevLex(vars) * DegRevLex(params)
    gb = groebner(sys, ordering=ordering)
    staircase = map(f -> Groebner.leading_term(f, ordering=ordering), gb)
    closed_staircase = all(var -> any(lead -> is_univariate_term(lead, var, vars), staircase), vars)
    params_relations = filter(poly -> all(x -> degree(poly, x) == 0, vars), gb)
    zerodim = closed_staircase && isempty(params_relations)
    zerodim
end

function nemo_crude_evaluate(poly::Nemo.FracElem, varmap)
    nemo_crude_evaluate(numerator(poly), varmap) // nemo_crude_evaluate(denominator(poly), varmap)
end

function nemo_crude_evaluate(poly::Rational{T}, varmap) where {T}
    poly
end

function nemo_crude_evaluate(poly::Nemo.FinFieldElem, varmap)
    poly
end

function nemo_crude_evaluate(poly::Nemo.ZZRingElem, varmap)
    BigInt(poly)
end

function nemo_crude_evaluate(poly::Nemo.MPolyRingElem, varmap)
    new_poly = 0
    for (i, term) in enumerate(Nemo.terms(poly))
        new_term = nemo_crude_evaluate(Nemo.coeff(poly, i), varmap)
        for var in Nemo.vars(term)
            exp = Nemo.degree(term, var)
            exp == 0 && continue
            new_var = varmap[var]
            new_term *= new_var^exp
        end
        new_poly += new_term
    end
    new_poly
end

function demote_gb(gb, nemo_vars, nemo_params)
    ring_flat = parent(nemo_vars[1])
    ring_param, params_demoted = Nemo.polynomial_ring(Nemo.base_ring(ring_flat), map(string, nemo_params))
    ring_demoted, vars_demoted = Nemo.polynomial_ring(Nemo.fraction_field(ring_param), map(string, nemo_vars), internal_ordering=:degrevlex)
    varmap = Dict((nemo_vars .=> vars_demoted)..., (nemo_params .=> params_demoted)...)
    gb_demoted = map(f -> ring_demoted(nemo_crude_evaluate(f, varmap)), gb)
    gb_result = empty(gb_demoted)
    without_i(arr, i) = arr[1:length(arr) .!= i]
    while true
        gb_demoted = sort(gb_demoted, by=leading_monomial)
        gb_demoted = map(f -> Nemo.map_coefficients(c -> c // Nemo.leading_coefficient(f), f), gb_demoted)
        gb_result = copy(gb_demoted)
        for i in 1:length(gb_demoted)
            f = gb_demoted[i]
            f_nf = Nemo.normal_form(f, without_i(gb_result, i))
            gb_result[i] = f_nf
            iszero(f_nf) && break
        end
        isequal(gb_demoted, gb_result) && break
        gb_result = filter(!iszero, gb_result)
        gb_demoted = gb_result
    end
    @assert all(f -> isone(Nemo.leading_coefficient(f)), gb_result)
    varmap_inv = Dict(v => k for (k, v) in varmap)
    gb_result = map(f -> nemo_crude_evaluate(f, varmap_inv), gb_result)
    gb_result = map(f -> f * denominator(f), gb_result)
    gb_result = gb_result .// one(ring_flat)
    @assert all(f -> total_degree(ring_flat(denominator(f))) == 0, gb_result)
    gb_result = map(numerator, gb_result)
    gb_result
end

function discriminant_variety_generically_zerodim(sys, vars, params)
    @assert internal_ordering(parent(sys[1])) == :degrevlex
    @assert is_generically_zerodim(sys, vars, params)
    
    # The only minor is the determinant
    minors = jacobian_minors(sys, vars, length(vars))
    W_c = [eliminate(vcat(sys, minors), vars)]

    ordering = DegRevLex(vars) * DegRevLex(params)
    gb = groebner(sys, ordering=ordering)
    gb = demote_gb(gb, vars, params) 
    
    staircase = map(f -> monomial(remove_vars(Groebner.leading_term(f, ordering=ordering), params), 1), gb)
    staircase_coeffs = []
    for (poly, lead) in zip(gb, staircase)
        cfs = filter(t -> lead == monomial(remove_vars(t, params), 1), collect(terms(poly)))
            cfs = map(c -> remove_vars(c, vars), cfs)
        push!(staircase_coeffs, sum(cfs))
    end
    W_infty = map(l -> [l], staircase_coeffs)

    W_d = vcat(W_c, W_infty)
    
    W_d
end

# Make it look more nice.
function postprocess(W_d; make_squarefree=true)
    W_d = filter(component -> !any(f -> is_constant(f), component), W_d)
    W_d_ = empty(W_d)
    for component in W_d
        if length(component) == 1
            factors = collect(factor(component[1]))
            for (factor, deg) in factors
                if make_squarefree
                    push!(W_d_, [factor])
                else
                    push!(W_d_, [factor^deg])
                end
            end
        else
            push!(W_d_, component)
        end
    end
    W_d = W_d_
    W_d = sort(W_d, by=f -> maximum(total_degree, f))
    W_d = unique(W_d)
    W_d
end

"""
    discriminant_variety(sys, vars, params)

Computes a Discriminant Variety of system `sys` from `Q[params][vars]`.
"""
function discriminant_variety(sys, vars, params)
    # Sanity checks
    @assert issubset(gens(parent(sys[1])), union(vars, params))
    @assert isempty(intersect(vars, params))

    ring_orig = parent(sys[1])
    K_orig = base_ring(ring_orig)
    @assert K_orig in (Nemo.QQ, AbstractAlgebra.QQ) || K_orig isa Nemo.FinField
    
    # Convert to Nemo in deg-rev-lex
    K = K_orig isa Nemo.FinField ? Nemo.GF(characteristic(K_orig)) : Nemo.QQ
    ring_drl, xs = polynomial_ring(K, symbols(ring_orig), internal_ordering=:degrevlex)
    varmap = Dict(gens(ring_orig) .=> gens(ring_drl))
    sys = map(f -> nemo_crude_evaluate(f, varmap), sys)
    vars = map(f -> nemo_crude_evaluate(f, varmap), vars)
    params = map(f -> nemo_crude_evaluate(f, varmap), params)

    if is_generically_zerodim(sys, vars, params)
        W_d = discriminant_variety_generically_zerodim(sys, vars, params)
    else
        error("Non-zerodim systems are not supported")
    end

    W_d = postprocess(W_d)
    
    W_d = map(F -> map(f -> change_base_ring(K_orig, f, parent=ring_orig), F), W_d)

    W_d
end

export discriminant_variety

### Example ###

#=
# Example 1 and Example 2 from
# https://members.loria.fr/GMoroz/assets/LGJMcca09.pdf

using Nemo

R, (x,y,z,a,b,c) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"])
sys = [a*x^2 + b - 1, y + b*z, y + c*z]
@show discriminant_variety(sys, [x,y,z], [a,b,c])

R, (x,y,a,b) = polynomial_ring(QQ, ["x", "y", "a", "b"])
sys = [a*x^6 + b*y^2 - 1, x^2 - a*y - b]
@show discriminant_variety(sys, [x,y], [a,b])
=#

end # module DiscriminantVariety
