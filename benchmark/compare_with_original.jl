# Compares `discriminant_variety` between a reference revision (default: HEAD,
# the original code) and the working tree (the new code).
#
# Usage, from anywhere:
#
#     julia benchmark/compare_with_original.jl [--tripod] [--timeout SECONDS] [--ref REV]
#
#   --tripod     also run Systems/tripod_1_1_qz.jl (the original code may need
#                a lot of time and memory on it)
#   --timeout    time limit per case and per version, in seconds (default 600)
#   --ref        git revision of the original code (default HEAD)
#
# Each version runs in its own temporary environment (fresh resolution of its
# Project.toml, the Manifest.toml of the repository is not used), and each case
# in its own Julia process. The time reported is the one of the call to
# discriminant_variety, after a warm-up on a small system.

const REPO = dirname(@__DIR__)

const CASES = [
    # name, file defining `sys` (or nothing), code building (sys, vars, params)
    ("basic", """
        R, (x, a) = polynomial_ring(QQ, ["x", "a"])
        ([a*x^2 + a + 1], [x], [a])"""),
    ("moroz_example_1", """
        R, (x, y, z, a, b, c) = polynomial_ring(QQ, ["x", "y", "z", "a", "b", "c"])
        ([a*x^2 + b - 1, y + b*z, y + c*z], [x, y, z], [a, b, c])"""),
    ("moroz_example_2", """
        R, (x, y, a, b) = polynomial_ring(QQ, ["x", "y", "a", "b"])
        ([a*x^6 + b*y^2 - 1, x^2 - a*y - b], [x, y], [a, b])"""),
    ("moroz_example_2_rational", """
        R, (x, y, a, b) = polynomial_ring(QQ, ["x", "y", "a", "b"])
        ([QQ(1,7)*a*x^6 + b*y^2 - QQ(1,3), x^2 - QQ(1,5)*a*y - b], [x, y], [a, b])"""),
    ("stability_bouzidi_rouillier", """
        include(joinpath(SYSTEMS, "stability_bouzidi_rouillier.jl"))
        xs = gens(parent(sys[1]))
        (sys, xs[1:2], xs[3:4])"""),
]
const TRIPOD = ("tripod_1_1_qz", """
        include(joinpath(SYSTEMS, "tripod_1_1_qz.jl"))
        xs = gens(parent(sys[1]))
        (sys, xs[1:5], xs[6:8])""")

# ---------------------------------------------------------------- worker
# julia --project=ENV compare_with_original.jl --worker SYSTEMS CASE_CODE OUTFILE
if length(ARGS) >= 1 && ARGS[1] == "--worker"
    systems, code, outfile = ARGS[2], ARGS[3], ARGS[4]
    @eval Main begin
        using Nemo, DiscriminantVariety
        const SYSTEMS = $systems
        # warm-up (compilation)
        let
            R, (x, a) = polynomial_ring(QQ, ["x", "a"])
            discriminant_variety([a*x^2 + a + 1], [x], [a])
        end
        sys_, vars_, params_ = eval(Meta.parse("begin\n" * $code * "\nend"))
        t = @elapsed W = discriminant_variety(sys_, vars_, params_)
        comps = sort([join(sort(map(string, C)), " ; ") for C in W])
        open($outfile, "w") do io
            println(io, t)
            foreach(c -> println(io, c), comps)
        end
    end
    exit(0)
end

# ---------------------------------------------------------------- driver
function parse_args(args)
    opts = Dict{String, Any}("tripod" => false, "timeout" => 600.0, "ref" => "HEAD")
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--tripod"
            opts["tripod"] = true
        elseif a == "--timeout"
            opts["timeout"] = parse(Float64, args[i += 1])
        elseif a == "--ref"
            opts["ref"] = args[i += 1]
        else
            error("unknown option $a")
        end
        i += 1
    end
    opts
end

julia() = Base.julia_cmd()

function make_env(dir, label)
    println("[$label] resolving and precompiling the environment in $dir ...")
    rm(joinpath(dir, "Manifest.toml"), force = true)
    run(`$(julia()) --project=$dir -e "import Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()"`)
end

function run_case(envdir, systems, code, timeout)
    out = tempname()
    log = tempname()
    cmd = `$(julia()) --project=$envdir $(@__FILE__) --worker $systems $code $out`
    p = run(pipeline(cmd, stdout = log, stderr = log), wait = false)
    t0 = time()
    while process_running(p) && time() - t0 < timeout
        sleep(0.2)
    end
    if process_running(p)
        kill(p)
        return (status = "TIMEOUT", time = NaN, result = String[])
    end
    if !success(p) || !isfile(out)
        lines = readlines(log)
        msg = isempty(lines) ? "failed" : last(filter(!isempty, lines), 3)
        return (status = "ERROR: " * join(msg, " | "), time = NaN, result = String[])
    end
    lines = readlines(out)
    (status = "ok", time = parse(Float64, lines[1]), result = lines[2:end])
end

function main(args)
    opts = parse_args(args)
    cases = opts["tripod"] ? vcat(CASES, [TRIPOD]) : CASES
    tmp = mktempdir()
    orig, new = joinpath(tmp, "original"), joinpath(tmp, "new")
    mkpath(orig)
    run(pipeline(`git -C $REPO archive --format=tar $(opts["ref"])`, `tar -x -C $orig`))
    mkpath(new)
    for f in ("Project.toml", "src", "Systems")
        cp(joinpath(REPO, f), joinpath(new, f))
    end
    make_env(orig, "original ($(opts["ref"]))")
    make_env(new, "new (working tree)")
    systems = joinpath(REPO, "Systems")
    println()
    rows = []
    for (name, code) in cases
        println("running $name ...")
        a = run_case(orig, systems, code, opts["timeout"])
        b = run_case(new, systems, code, opts["timeout"])
        same = a.status == "ok" && b.status == "ok" ? (a.result == b.result ? "yes" : "NO") : "-"
        push!(rows, (name, a, b, same))
    end
    fmt(r) = r.status == "ok" ? string(round(r.time, digits = 3), " s") : r.status
    println("\n", rpad("case", 30), rpad("original", 22), rpad("new", 22), "same result")
    for (name, a, b, same) in rows
        println(rpad(name, 30), rpad(fmt(a), 22), rpad(fmt(b), 22), same)
    end
    for (name, a, b, same) in rows
        same == "NO" || continue
        println("\n$name: components differ")
        println("  original: ", join(a.result, "\n            "))
        println("  new:      ", join(b.result, "\n            "))
    end
    rm(tmp, recursive = true, force = true)
end

main(ARGS)
