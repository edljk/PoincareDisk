module PoincareDisk 

using LinearAlgebra, StatsBase, Random, CoordinateTransformations, Roots, Printf
using Makie.GeometryBasics, ColorTypes

include("utils.jl")

"""
   C, r, P, Q, arclength = geodesic𝔻(U::Array{T}, V::Array{T}) where T

Find the geodesic arc through the two points U and V perpendicular to the unit circle

N.B. See 
https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model#Analytic_geometry_constructions_in_the_hyperbolic_plane
"""
function geodesic𝔻(U::Array{T}, V::Array{T}) where T
    ε, detUV = T(10) ^ (- 10), U[1] * V[2] - U[2] * V[1]
    n2U, n2V = U[1] ^ 2 + U[2] ^ 2, V[1] ^ 2 + V[2] ^ 2
    if n2U ≥ 1 || n2V ≥ 1 
        error("pb of data points in geodesic (colinear or not in 𝔻)")
    end
    # colinear configuration
    if abs(detUV) < ε
        C, r = zeros(T, 2), T(Inf)
        P = (U - V) / norm(U - V)
        Q = - P
        arclength = log((1 + norm(U)) / (1 - norm(U))) + 
                    log((1 + norm(V)) / (1 - norm(V)))
        return C, r, P, Q, arclength
    end
    a = (U[2] * (n2V + 1) - V[2] * (n2U + 1)) / detUV
    b = (V[1] * (n2U + 1) - U[1] * (n2V + 1)) / detUV
    C = - [a / 2, b / 2]
    r = sqrt(- 1 + (a ^ 2 + b ^ 2) / 4)
    # intersection between the two circles
    P, Q = intersect_circles(zeros(T, 2), T(1), C, r) 
    # geodesic distance 
    fP = max(norm(U - P) / norm(V - P), norm(V - P) / norm(U - P)) 
    fQ = max(norm(U - Q) / norm(V - Q), norm(V - Q) / norm(U - Q))  
    arclength = log(fP * fQ)
    return C, r, P, Q, arclength
end
#-------------------------------------------------------------------------------
""" 
   C, R = circle𝔻(I::Array{T}, r::T) where T

Associate to every radius / center in 𝔻 its euclidean representation
"""
function circle𝔻(I::Array{T}, r::T) where T
    nI = norm(I)
    d = I / nI
    dOI = 2 * atanh(nI)
    #_, _, _, _, dOI = geodesic𝔻(zeros(2), I)
    OA, OB = tanh((dOI + r) / 2), tanh((dOI - r) / 2)
    R = (OA - OB) / 2
    OIp = (OA + OB) / 2
    return OIp * d, R
end
#-------------------------------------------------------------------------------
""" 
   P, Q = equitriangle𝔻(U::Array{T}, V::Array{T}) where T

Identify the two equilaterial triangles associated to [U, V] ∈ 𝔻
"""
function equitriangle𝔻(U::Array{T}, V::Array{T}) where T
    _, _, _, _, dUV = geodesic𝔻(U, V)
    CU, rU = circle𝔻(U, dUV)   
    CV, rV = circle𝔻(V, dUV)   
    P, Q = intersect_circles(CU, rU, CV, rV) 
    return P, Q
end
#-------------------------------------------------------------------------------
"""
   C, r, Fd = perpendicular_by_P𝔻(U::Array{T}, V::Array{T}, M::Array{T}, d::T) 

Compute the geodesic perpendicular arc passing throuch M ∈ [U,V]. Fd lies on 
this arc at distance d of M
"""
function perpendicular_by_P𝔻(U::Array{T}, V::Array{T}, M::Array{T}, 
                             d::T = 0) where T 
    ε = 1e-12
    if abs(U[1] * V[2] - U[2] * V[1]) < ε # almost aligned points
        a = norm(M)
        b = a + (1 - a ^ 2) / (2 * a)
        Cper, rper = b * M / norm(M), b - a
        CM, rM = circle𝔻(M, d)    
        Fd = intersect_circles(CM, rM, Cper, rper)[1]
        return Cper, rper, Fd 
    end
    C, r, P, Q, lUV = geodesic𝔻(U, V)

    lMU = norm(M - U) > ε ? geodesic𝔻(M, U)[5] : 0.
    lMV = norm(M - V) > ε ? geodesic𝔻(M, V)[5] : 0.
    rε = abs(d) < 1e-12 ? min(lMU, lMV) / 4. : d 
    CM, rM = circle𝔻(M, rε)
    A, B = intersect_circles(C, r, CM, rM)
    C1, C2 = equitriangle𝔻(A, B)
    F = norm(C1 - C) < norm(C2 - C) ? C2 : C1
    Cper, rper, _ = geodesic𝔻(F, M)
    E, F = intersect_circles(Cper, rper, CM, rM)
    Fd = norm(E - C) < r ? F : E
    return Cper, rper, Fd
end
#-------------------------------------------------------------------------------
"""
   Plist = hexagon𝔻(l::Vector{T}, yinit::T = T(5))

Find (centered) vertices of an hexagon of prescribed alternate lengths
"""
function hexagon𝔻(l::Vector{T}, yinit::T = T(10)) where T
    # start with a large value
    y = yinit
    f(x) = hexagon_cost(x, l)[2] - l[3]
    # search root 
    atol = typeof(yinit) == BigFloat ? 1e-20 : 1e-8
    rtol = atol
    sols = find_zeros(f, (yinit / 100, yinit), atol = atol, rtol = rtol)
    sol = sols[findmin(abs.(f.(sols)))[2]]
    # reconstruct hexagon 
    Plist, _ = hexagon_cost(sol, l, debug = false)
    return Plist
end
"""
    Plist = hexagon_cost(y::T) where T

"""
function hexagon_cost(y::T, l::Vector{T}; 
                      debug::Bool = false) where T
    Plist = Vector{Vector{T}}()
    P1 = zeros(T, 2)
    push!(Plist, P1)
    P2 = [(exp(l[1]) - 1) / (exp(l[1]) + 1), T(0)]
    push!(Plist, P2)
    P3 = perpendicular_by_P𝔻(P1, P2, P2, y)[3]
    # stay in the non negative cadran
    if P3[2] < 0
        P3[2] *= - 1
    end
    push!(Plist, P3)
    P4 = perpendicular_by_P𝔻(P2, P3, P3, l[2])[3]
    push!(Plist, P4)
    C4, r4, _ = perpendicular_by_P𝔻(P3, P4, P4, T(0.1))
    # the disk must fully be contained in x > 0
    if minimum(C4 .- r4) < 0
        if debug 
            @warn "last computed geodesic not in x > 0"
        end
        return Plist, - Inf
    end
    # orthogonal geodesic WARNING: this step is only relevant in the specific 
    # context of the first ortho quarter
    b2 = C4[2]
    c2 = (r4 ^ 2 - norm(C4) ^ 2) / 2
    Corth = [0., (0.5 - c2) / b2]
    Rorth = sqrt(norm(Corth - C4) ^ 2 - r4 ^ 2)
    # last vertices
    A, B = intersect_circles(Corth, Rorth, C4, r4)
    I1 = norm(A) > 1 ? B : A
    push!(Plist, I1)
    D, E = intersect_linecircle([T(0), T(0)], [T(0), T(1)], Corth, Rorth)
    I2 = norm(D) > 1 ? E : D
    push!(Plist, I2)
    # 3rd length 
    l3 = geodesic𝔻(I1, I2)[5]
    return Plist, l3
end
end # module 