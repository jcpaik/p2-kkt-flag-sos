#!/usr/bin/env julia
# Frobenius-canonical certificate selector in arbitrary precision.
#
#   min 1/2 ||Y||_F^2   s.t.  F_i . Y = c_i  (i = 1..m),
#                             F_0 . Y >= t0,
#                             Y block-diagonal PSD,
#
# where (F_0, F_i, c) come from an exported PROBLEM.dat-s (the SDPA-P
# data whose dual side carries the certificate).  Unlike the trace
# selector, the objective is strictly convex, so the minimizer is unique
# and varies continuously with the data — the robust canonical point for
# cross-degree/eps regression.
#
# Usage: julia clarabel_selector.jl PROBLEM.dat-s t0 output.json

using Clarabel, SparseArrays, LinearAlgebra, GenericLinearAlgebra, JSON
using Printf

setprecision(BigFloat, 256)
T = BigFloat

function parse_dats(path::String)
    lines = readlines(path)
    m = parse(Int, split(lines[1], "=")[1])
    nblock = parse(Int, split(lines[2], "=")[1])
    struct_str = strip(split(lines[3], "=")[1])
    sizes = [parse(Int, strip(t)) for t in split(strip(struct_str, ['(', ')', ' ']), ",")]
    cline = strip(lines[4])
    cvals = [parse(T, strip(t)) for t in split(strip(cline, ['{', '}']), ",")]
    entries = Vector{Tuple{Int,Int,Int,Int,T}}()
    for line in lines[5:end]
        isempty(strip(line)) && continue
        parts = split(strip(line))
        push!(entries, (parse(Int, parts[1]), parse(Int, parts[2]),
                        parse(Int, parts[3]), parse(Int, parts[4]),
                        parse(T, parts[5])))
    end
    return m, sizes, cvals, entries
end

function main()
    problem, t0str, output = ARGS[1], ARGS[2], ARGS[3]
    t0 = parse(T, t0str)
    m, sizes, c, entries = parse_dats(problem)
    @printf("m=%d blocks=%d\n", m, length(sizes))

    # svec layout: blocks in order, column-major upper triangle per
    # Clarabel's PSDTriangleCone convention, off-diagonals scaled sqrt(2).
    offsets = Int[]
    total = 0
    for n in sizes
        push!(offsets, total)
        total += n * (n + 1) ÷ 2
    end
    svec_index(b, i, j) = begin  # requires i <= j (upper triangle, col-major)
        offsets[b] + (j - 1) * j ÷ 2 + i
    end
    rt2 = sqrt(T(2))

    # Constraint matrices: rows 1..m equalities (zero cone),
    # row m+1: -F0.Y <= -t0 (nonnegative cone).
    rowsE = Int[]; colsE = Int[]; valsE = T[]
    rows0 = Int[]; cols0 = Int[]; vals0 = T[]
    for (k, b, i, j, v) in entries
        i2, j2 = min(i, j), max(i, j)
        col = svec_index(b, i2, j2)
        w = i2 == j2 ? v : rt2 * v
        if k == 0
            push!(rows0, 1); push!(cols0, col); push!(vals0, w)
        else
            push!(rowsE, k); push!(colsE, col); push!(valsE, w)
        end
    end

    A_eq = sparse(rowsE, colsE, valsE, m, total)
    # F0 row as a dense-ish sparse row; constraint -F0.Y + s = -t0, s>=0.
    A_f0 = sparse(rows0, cols0, -vals0, 1, total)
    # PSD cones: -Y + s = 0, s in PSD triangle.
    A_psd = -sparse(I, total, total) * one(T)

    A = [A_eq; A_f0; A_psd]
    b = [c; -t0; zeros(T, total)]
    cones = Clarabel.SupportedCone[Clarabel.ZeroConeT(m),
                                    Clarabel.NonnegativeConeT(1)]
    for n in sizes
        push!(cones, Clarabel.PSDTriangleConeT(n))
    end

    P = sparse(I, total, total) * one(T)
    q = zeros(T, total)

    settings = Clarabel.Settings{T}(
        verbose = true,
        max_iter = 200,
        tol_gap_abs = T(1e-25),
        tol_gap_rel = T(1e-25),
        tol_feas = T(1e-25),
    )
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, A, b, cones, settings)
    solution = Clarabel.solve!(solver)

    x = solution.x
    frob2 = dot(x, x)
    # Trace: sum of diagonal svec entries.
    tracev = zero(T)
    for (bidx, n) in enumerate(sizes)
        for d in 1:n
            tracev += x[svec_index(bidx, d, d)]
        end
    end
    @printf("status=%s  ||Y||_F^2=%.10e  trace=%.10e\n",
            string(solution.status), Float64(frob2), Float64(tracev))

    blocks = Vector{Vector{Vector{String}}}()
    for (bidx, n) in enumerate(sizes)
        M = Matrix{T}(undef, n, n)
        for j in 1:n, i in 1:j
            v = x[svec_index(bidx, i, j)]
            v = i == j ? v : v / rt2
            M[i, j] = v; M[j, i] = v
        end
        push!(blocks, [[string(M[i, j]) for j in 1:n] for i in 1:n])
    end
    open(output, "w") do io
        JSON.print(io, Dict(
            "status" => string(solution.status),
            "frobenius_sq" => string(frob2),
            "trace" => string(tracev),
            "t0" => t0str,
            "blocks" => blocks,
        ))
    end
    println("wrote ", output)
end

main()
