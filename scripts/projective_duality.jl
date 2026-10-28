"""
Projective-dual degrees via Oscar's IntersectionTheory.

Tevelev, Theorem 7.2 (arXiv version): put
    f(q) = q integral(c_top(Omega_X tensor O_X(q))).
The first nonconstant term at q=1 has exponent codim(X^dual), and its
coefficient is deg(X^dual). Equivalently, Theorem 7.1 gives the jet-bundle
Chern numbers recorded below. The ample generator is the fundamental
weight omega_k, including in the spinor cases.

Run: julia --project=scripts scripts/projective_duality.jl [--check]
     julia --project=scripts scripts/projective_duality.jl --case C6/3
"""

using Oscar
import JSON

const ROOT = dirname(@__DIR__)
const TARGET = joinpath(ROOT, "references", "projective-duality.json")

function geometry(T, n, k)
    # Tevelev, Theorem 9.19: the defective maximal-parabolic embeddings.
    empty = (T == :A && k in (1, n)) || (T == :C && k == 1) ||
            (T, n, k) == (:B, 2, 2) || (T == :D && n == 3 && k in (2, 3))
    defect = T == :A && n >= 4 && iseven(n) && k in (2, n - 1) ? 2 : 0
    (T, n, k) in ((:B, 4, 4), (:D, 5, 4), (:D, 5, 5)) && (defect = 4)
    return empty, defect
end

coordinates(x, n) = ntuple(i -> Int(coefficients(x)[1, i]), n)

function root_data(T, n, k)
    # Oscar reserves type D for rank >= 4; retain the catalogue's D3
    # numbering (node 1 is central, nodes 2 and 3 are the spinor nodes).
    R = T == :D && n == 3 ? root_system([2 -1 -1; -1 2 0; -1 0 2]) : root_system(T, n)
    roots = positive_roots(R)
    root_coordinates = [coordinates(alpha, n) for alpha in roots]
    coroots = [coordinates(dual(alpha), n) for alpha in roots]
    simple = [coordinates(WeightLatticeElem(alpha), n) for alpha in simple_roots(R)]
    weights = [coordinates(mu, n) for mu in weyl_orbit(fundamental_weight(R, k))]
    d = count(alpha -> alpha[k] != 0, root_coordinates)
    return (; n, roots=root_coordinates, coroots, simple, weights, d)
end

dot_product(a, b) = sum(a[i] * b[i] for i in eachindex(a))

function chern_numbers(data, defect, parameters)
    (; n, roots, coroots, simple, weights, d) = data
    # A regular integral cocharacter, specified by its simple-root values.
    A = matrix(QQ, n, n, [simple[i][j] for i in 1:n for j in 1:n])
    point = inv(A) * matrix(QQ, n, 1, parameters)
    scale = lcm([denominator(point[i, 1]) for i in 1:n])
    values = [ZZ(scale * point[i, 1]) for i in 1:n]
    simple_values = [dot_product(alpha, values) for alpha in simple]
    root_values = [dot_product(alpha, simple_values) for alpha in roots]
    @assert all(>(0), root_values)

    # At mu=w(omega_k), cotangent roots are alpha with <mu,alpha^vee> > 0.
    # The fiber of the ample generator has weight -mu.
    function cotangent_weights(mu)
        result = ZZRingElem[]
        for (alpha, coroot) in zip(root_values, coroots)
            pairing = dot_product(mu, coroot)
            pairing == 0 || push!(result, pairing > 0 ? alpha : -alpha)
        end
        @assert length(result) == d
        return result
    end

    X = tn_variety(d, [mu => 1 for mu in weights])
    X.T = tn_bundle(X, d, mu -> TnRep(-cotangent_weights(mu)))
    L = tn_bundle(X, 1, mu -> TnRep([-dot_product(mu, values)]))
    degree_X = integral(chern_class(L, 1)^d)

    # Taylor-expand at q=1. Truncation keeps even E8/P4 (dimension 106)
    # practical, while Oscar performs both the Chern-class and integral steps.
    S, t = power_series_ring(QQ, defect + 2, "t")
    q = 1 + t
    twisted_cotangent = tn_bundle(X, d, mu -> begin
        shift = q * dot_product(mu, values)
        TnRep([S(alpha) - shift for alpha in cotangent_weights(mu)])
    end)
    f = q * integral(top_chern_class(twisted_cotangent))
    numbers = [degree_X; [coeff(f, j) for j in 1:defect + 1]]
    @assert all(is_integer, numbers)
    return BigInt.(ZZ.(numbers))
end

function compute(row)
    label, k = row["type"], row["parabolic"]
    T, n = Symbol(label[1:1]), parse(Int, label[2:end])
    data = root_data(T, n, k)
    empty, defect = geometry(T, n, k)
    first = chern_numbers(data, defect, ones(Int, n))
    second = chern_numbers(data, defect, collect(1:n))
    @assert first == second "cocharacter dependence: $label/$k"
    @assert first[1] == row["degree"] "incorrect polarization: $label/$k"
    @assert all(iszero, first[2:end-1]) "incorrect dual defect: $label/$k"
    @assert empty ? iszero(first[end]) : first[end] > 0
    return Dict(
        "dimension" => data.d, "fixed_points" => length(data.weights),
        "degree_of_grassmannian" => first[1],
        "codimension" => empty ? nothing : defect + 1,
        "degree" => empty ? nothing : first[end], "jet_integrals" => first[2:end],
    )
end

function main(args)
    @assert args in (String[], ["--check"]) || (length(args) == 2 && args[1] == "--case") "Use --check or --case C6/3"
    selected = length(args) == 2 ? args[2] : nothing
    rows = JSON.parsefile(joinpath(ROOT, "data", "grassmannians.json"))
    results = Dict{String, Any}()
    for row in rows
        parse(Int, row["type"][2:end]) > 8 && continue
        key = "$(row["type"])/$(row["parabolic"])"
        selected !== nothing && key != selected && continue
        started = time()
        results[key] = compute(row)
        println(key, " ", results[key], " (", round(time() - started; digits=2), "s)")
        flush(stdout)
    end
    @assert !isempty(results) "No catalogue entry: $selected"
    if selected !== nothing
        return
    elseif args == ["--check"]
        @assert results == JSON.parsefile(TARGET)
        println("All $(length(results)) records agree with both exact Oscar evaluations.")
    else
        open(TARGET, "w") do io
            JSON.print(io, results, 2)
            println(io)
        end
    end
end

main(ARGS)
