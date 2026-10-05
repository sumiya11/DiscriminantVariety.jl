# Multi-modular computations for the block ordering DRL(vars) > DRL(params),
# lifting over QQ only what the discriminant variety needs.
#
# Notation: I = <sys> in QQ[vars, params].
#
#   W_c, W_s, saturations: the elimination ideal I ∩ QQ[params];
#   W_inf:                 for each variable x of `vars`, the leading coefficients
#                          (polynomials in `params`) of the elements of a Groebner
#                          basis for DRL(vars) > DRL(params) whose leading monomial
#                          in `vars` is a pure power of x. These coefficients form a
#                          (non reduced) Groebner basis in the parameters.
#
# Neither needs the reduced Groebner basis of I, which is large and expensive to
# obtain (most of the time is then spent in the final reduction). The method:
#
#   1. Homogenize the system with a new variable h, placed last in the block of the
#      parameters: DRL(vars) > DRL(params, h). No saturation by h is needed: if G is
#      a Groebner basis of the homogenized system, G with h = 1 is a (non reduced)
#      Groebner basis of I for DRL(vars) > DRL(params). On homogeneous input the
#      ordinary degree selection of F4 is the right one for a block ordering.
#   2. Modulo a prime p: Groebner basis of the homogenized system, recording its
#      trace (learn, trace 1). Keep only the part of interest, with h = 1.
#   3. Reduced Groebner basis of each kept ideal, a small ideal in the parameters
#      (learn, trace 2).
#   4. For further primes: apply trace 1, keep, apply trace 2. The coefficients of
#      the kept part only are lifted to QQ (CRT + rational reconstruction), until
#      a new prime confirms the result.
#
# Only the documented low-level interface of Groebner.jl is used (exponent vectors
# and coefficients modulo p): groebner_learn, groebner_apply!, Groebner.PolyRing.

const _G = Groebner
const _NZ = Nemo.ZZ

# The input in the format of the low-level interface of Groebner.jl: exponent
# vectors in the variables (vars..., params..., h), homogenized, and the
# coefficients as numerators and denominators.
struct HomogenizedBlockSystem
    nx::Int                                  # number of eliminated variables
    nu::Int                                  # number of parameters
    exps::Vector{Vector{Vector{Int}}}
    num::Vector{Vector{BigInt}}
    den::Vector{Vector{BigInt}}
end

function HomogenizedBlockSystem(sys, vars, params)
    R = parent(first(sys))
    position = Dict(g => i for (i, g) in enumerate(gens(R)))
    perm = vcat([position[v] for v in vars], [position[u] for u in params])
    exps, num, den = Vector{Vector{Vector{Int}}}(), Vector{Vector{BigInt}}(), Vector{Vector{BigInt}}()
    for f in sys
        iszero(f) && continue
        E = [Int.(e[perm]) for e in exponent_vectors(f)]
        d = maximum(sum, E)
        push!(exps, [vcat(e, d - sum(e)) for e in E])
        cs = collect(coefficients(f))
        push!(num, [BigInt(numerator(c)) for c in cs])
        push!(den, [BigInt(denominator(c)) for c in cs])
    end
    HomogenizedBlockSystem(length(vars), length(params), exps, num, den)
end

# Coefficients modulo p, or nothing if p divides a numerator or a denominator.
function residues(S::HomogenizedBlockSystem, p)
    out = Vector{Vector{UInt64}}(undef, length(S.num))
    for i in eachindex(S.num)
        v = Vector{UInt64}(undef, length(S.num[i]))
        for j in eachindex(v)
            n, d = mod(S.num[i][j], p), mod(S.den[i][j], p)
            (iszero(n) || iszero(d)) && return nothing
            v[j] = UInt64(mod(n * invmod(d, p), p))
        end
        out[i] = v
    end
    out
end

block_ring(S::HomogenizedBlockSystem, p) =
    _G.PolyRing(S.nx + S.nu + 1, DegRevLex(collect(1:S.nx)) * DegRevLex(collect(S.nx+1:S.nx+S.nu+1)), p)
param_ring(S::HomogenizedBlockSystem, p) = _G.PolyRing(S.nu, DegRevLex(), p)

# The part of interest of a basis (gbm, gbc) of the homogenized system, with h = 1:
# a list of ideals in the parameters, each as (exponent vectors, coefficients).
#   :elim -> one ideal: the elements free of the eliminated variables;
#   :winf -> one ideal per eliminated variable: the leading coefficients of the
#            elements whose leading monomial in `vars` is a pure power of it.
# The terms of each element are sorted by decreasing monomial: the first one is
# the leading term.
function kept_part(mode, S::HomogenizedBlockSystem, gbm, gbc)
    nx, nv = S.nx, S.nx + S.nu + 1
    par(e) = e[nx+1:nv-1]
    xpart(e) = view(e, 1:nx)
    if mode === :elim
        idx = [i for i in eachindex(gbm) if all(e -> all(iszero, xpart(e)), gbm[i])]
        return [([map(par, gbm[i]) for i in idx], [gbc[i] for i in idx])]
    end
    @assert mode === :winf
    out = []
    for v in 1:nx
        M, C = Vector{Vector{Vector{Int}}}(), Vector{Vector{eltype(gbc[1])}}()
        for i in eachindex(gbm)
            m = xpart(gbm[i][1])
            (m[v] > 0 && all(w -> w == v || iszero(m[w]), 1:nx)) || continue
            sel = [j for j in eachindex(gbm[i]) if xpart(gbm[i][j]) == m]
            push!(M, [par(gbm[i][j]) for j in sel])
            push!(C, gbc[i][sel])
        end
        push!(out, (M, C))
    end
    out
end

# Is the system generically zero-dimensional in `vars`, read on the basis modulo p:
# every eliminated variable has a pure power as leading monomial, and no element is
# free of the eliminated variables.
function zerodim_from_basis(S::HomogenizedBlockSystem, gbm)
    nx = S.nx
    lead(i) = view(gbm[i][1], 1:nx)
    closed = all(v -> any(i -> lead(i)[v] > 0 && all(w -> w == v || iszero(lead(i)[w]), 1:nx), eachindex(gbm)), 1:nx)
    relations = any(i -> all(e -> all(iszero, view(e, 1:nx)), gbm[i]), eachindex(gbm))
    closed && !relations
end

is_unit_ideal(M) = any(m -> length(m) == 1 && all(iszero, m[1]), M)

# Primes below 2^31, by decreasing value.
mutable struct PrimeSequence
    p::Int
end
PrimeSequence() = PrimeSequence(2^31 + 1)
function next_prime!(s::PrimeSequence)
    p = s.p - 2
    while !Nemo.is_prime(p)
        p -= 2
    end
    s.p = p
end

# One prime for which the coefficients of the system are invertible.
function next_good_prime!(s::PrimeSequence, S)
    while true
        p = next_prime!(s)
        r = residues(S, p)
        r === nothing || return p, r
    end
end

# Sizes of batches accepted by groebner_apply!.
const BATCH_SIZES = (1, 2, 4, 8, 16, 32, 64, 128)

# Applies trace 1 modulo the primes of `batch` (pairs (p, residues)), together when
# there are several of them. Returns the pairs (p, coefficients) of the primes that
# follow the trace. A batch in which one prime fails is replayed prime by prime, to
# drop only that prime.
function apply_trace1(tr1, S::HomogenizedBlockSystem, batch)
    if length(batch) > 1
        ok, cs = groebner_apply!(tr1, Tuple((block_ring(S, p), S.exps, r) for (p, r) in batch))
        ok && return [(batch[k][1], cs[k]) for k in eachindex(batch)]
    end
    out = []
    for (p, r) in batch
        ok, c = groebner_apply!(tr1, block_ring(S, p), S.exps, r)
        ok && push!(out, (p, c))
    end
    out
end

"""
    lift_kept_part(S, mode; zerodim_check = false, batch = 4, verbose = false, maxprimes = 2000)

Steps 2 to 4 of the method (see the top of this file). Returns `(zerodim, ideals)`:
`ideals` is a list of ideals in the parameters, each as `(exponents, coefficients)`
with rational coefficients (its reduced Groebner basis), or `:unit` for the unit
ideal. With `zerodim_check = true`, the generic zero-dimensionality is first read on
the basis modulo the first prime; a negative answer is confirmed with a second
prime, and nothing is lifted if it is confirmed.

After the first prime, trace 1 is replayed for `batch` primes at a time (one of
$(BATCH_SIZES)); `batch = 1` replays one prime after the other. Within a batch the
primes are used one by one, and the lifting stops as soon as one confirms the result.
"""
function lift_kept_part(S::HomogenizedBlockSystem, mode; zerodim_check = false, batch = 4, verbose = false, maxprimes = 2000)
    batch in BATCH_SIZES || throw(ArgumentError("batch must be one of $(BATCH_SIZES)"))
    primes = PrimeSequence()
    p, res = next_good_prime!(primes, S)
    tr1, gbm, gbc = groebner_learn(block_ring(S, p), S.exps, res; homogenize = :no)
    if zerodim_check && !zerodim_from_basis(S, gbm)
        # A negative answer is confirmed with a second prime. If the second prime
        # disagrees, the first one was unlucky and the second one is used instead.
        p, res = next_good_prime!(primes, S)
        tr1, gbm, gbc = groebner_learn(block_ring(S, p), S.exps, res; homogenize = :no)
        zerodim_from_basis(S, gbm) || return false, nothing
    end
    K = kept_part(mode, S, gbm, gbc)
    ideals = Vector{Any}(undef, length(K))
    todo = Int[]
    for k in eachindex(K)
        if isempty(K[k][1])
            ideals[k] = (Vector{Vector{Int}}[], Vector{Rational{BigInt}}[])
        elseif is_unit_ideal(K[k][1])
            ideals[k] = :unit
        else
            push!(todo, k)
        end
    end
    isempty(todo) && return true, ideals
    # trace 2: reduced basis of each kept ideal
    tr2, sup, acc = Dict{Int, Any}(), Dict{Int, Any}(), Dict{Int, Vector{Vector{Nemo.ZZRingElem}}}()
    for k in todo
        tr2[k], sup[k], c = groebner_learn(param_ring(S, p), K[k][1], [UInt64.(x) for x in K[k][2]])
        acc[k] = [_NZ.(x) for x in c]
    end
    modulus = _NZ(p)
    candidate = nothing
    nprimes, nused = 1, 1
    while nprimes < maxprimes
        primes_batch = [next_good_prime!(primes, S) for _ in 1:batch]
        nprimes += batch
        for (p, c1) in apply_trace1(tr1, S, primes_batch)
            Kp = kept_part(mode, S, gbm, c1)
            cur = Dict{Int, Any}()
            for k in todo
                ok, c2 = groebner_apply!(tr2[k], param_ring(S, p), K[k][1], [UInt64.(x) for x in Kp[k][2]])
                ok || break
                cur[k] = c2
            end
            length(cur) == length(todo) || continue
            # the candidate is accepted when this new prime confirms it
            if candidate !== nothing && all(k -> agrees(candidate[k], cur[k], p), todo)
                verbose && println("lifted with $nused primes, checked with one more ($nprimes primes computed)")
                for k in todo
                    ideals[k] = (sup[k], candidate[k])
                end
                return true, ideals
            end
            for k in todo, i in eachindex(acc[k]), j in eachindex(acc[k][i])
                acc[k][i][j] = Nemo.crt(acc[k][i][j], modulus, _NZ(cur[k][i][j]), _NZ(p))
            end
            modulus *= p
            nused += 1
            candidate = reconstruct_all(acc, todo, modulus)
        end
    end
    error("the multi-modular lifting did not converge with $maxprimes primes")
end

function reconstruct_all(acc, todo, modulus)
    bound = isqrt(div(modulus, 2) - 1)
    out = Dict{Int, Vector{Vector{Rational{BigInt}}}}()
    for k in todo
        out[k] = Vector{Vector{Rational{BigInt}}}(undef, length(acc[k]))
        for i in eachindex(acc[k])
            out[k][i] = Vector{Rational{BigInt}}(undef, length(acc[k][i]))
            for j in eachindex(acc[k][i])
                ok, q = Nemo.reconstruct(acc[k][i][j], modulus, bound, bound)
                ok || return nothing
                out[k][i][j] = Rational{BigInt}(BigInt(numerator(q)), BigInt(denominator(q)))
            end
        end
    end
    out
end

agrees(cand, cur, p) = all(zip(cand, cur)) do (q, c)
    all(j -> mod(numerator(q[j]) * invmod(denominator(q[j]), p), p) == c[j], eachindex(q))
end

# Polynomials of `R` in the variables `params` from exponent vectors in `params`.
function to_polynomials(R, params, monoms, coeffs)
    position = Dict(g => i for (i, g) in enumerate(gens(R)))
    pos = [position[u] for u in params]
    K = base_ring(R)
    map(zip(monoms, coeffs)) do (M, C)
        ctx = MPolyBuildCtx(R)
        for (m, c) in zip(M, C)
            e = zeros(Int, ngens(R))
            e[pos] .= m
            push_term!(ctx, K(numerator(c)) // K(denominator(c)), e)
        end
        finish(ctx)
    end
end

"""
    eliminate_modular(sys, vars; batch = 4, verbose = false)

Reduced Groebner basis of the elimination ideal `<sys> ∩ QQ[params]`, where the
parameters are the other variables of the ring (see the top of this file).
`batch`: number of primes replayed together (`1`: one prime after the other).
"""
function eliminate_modular(sys, vars; batch = 4, verbose = false)
    R = parent(first(sys))
    params = filter(g -> !(g in vars), gens(R))
    isempty(params) && error("eliminate_modular: no variable left")
    S = HomogenizedBlockSystem(sys, vars, params)
    _, ideals = lift_kept_part(S, :elim; batch = batch, verbose = verbose)
    I = only(ideals)
    I === :unit && return [one(R)]
    to_polynomials(R, params, I[1], I[2])
end

"""
    infinity_modular(sys, vars, params; batch = 4, verbose = false) -> (zerodim, W_inf)

Checks that the system is generically zero-dimensional in `vars` (on a basis
modulo a prime) and, if so, returns the components of W_inf: for each eliminated
variable whose ideal of leading coefficients is not the unit ideal, the reduced
Groebner basis of this ideal (see the top of this file). `batch`: as for
`eliminate_modular`.
"""
function infinity_modular(sys, vars, params; batch = 4, verbose = false)
    R = parent(first(sys))
    S = HomogenizedBlockSystem(sys, vars, params)
    zerodim, ideals = lift_kept_part(S, :winf; zerodim_check = true, batch = batch, verbose = verbose)
    zerodim || return false, nothing
    W = Vector{Vector{elem_type(R)}}()
    for I in ideals
        (I === :unit || isempty(I[1])) && continue
        push!(W, to_polynomials(R, params, I[1], I[2]))
    end
    true, W
end
