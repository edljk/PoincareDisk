"""
   p = arcparam𝔻(U::Array{T}, V::Array{T}, C::Array{T}, r::T, np) where T

Discretization of the geodesic arc (UV) ∈ 𝔻
"""
function arcparam𝔻(U::Array{T}, V::Array{T}, C::Array{T}, r::T, np) where T
    if r == T(Inf) # diameter 
        t = collect(LinRange(0, 1, np))
        return t .* U' + (1 .- t ) .* V'
    end
    θU, θV = PolarFromCartesian()(U - C).θ, PolarFromCartesian()(V - C).θ
    # angles in [0, 2π]
    θU < 0 && (θU = 2 * T(π) + θU)
    θV < 0 && (θV = 2 * T(π) + θV)
    if θU > θV
        θU, θV = θV, θU 
    end
    l1, l2 = θV - θU, θU + 2 * T(π) - θV 
    # consider the arc of smallest arc 
    θ = if l1 < l2 
        LinRange(θU, θV, np)
    else
        LinRange(θV - 2 * T(π), θU, np)
    end
    p = C' .+ r * hcat(cos.(θ), sin.(θ))
    return p
end
#-------------------------------------------------------------------------------
"""
   
   P, Q = intersect_circles(C1::Array{T}, r1::T, C2::Array{T}, r2::T) where T

Intersection of the two euclidan disks (C1, r1) and (C2, r2)
"""
function intersect_circles(C1::Array{T}, r1::T, C2::Array{T}, r2::T) where T
    C = C2 - C1
    M, R = (C1 + C2) / 2, norm(C)
    if R > (r1 + r2)
        error("empty intersection in intersect_circles")
    end
    sqri = sqrt(2 * (r1^ 2 + r2 ^ 2) / R ^ 2 - (r1 ^ 2 - r2 ^ 2) ^ 2 / R ^ 4 -
                1) 
    P = M + (r1 ^ 2 - r2 ^ 2) / (2 * R ^ 2) * C + 0.5 * sqri* [C[2], - C[1]]
    Q = M + (r1 ^ 2 - r2 ^ 2) / (2 * R ^ 2) * C - 0.5 * sqri* [C[2], - C[1]]
    return P, Q
end
#-------------------------------------------------------------------------------
"""
    P, Q = intersect_linecircle(P1::Array{T}, P2::Array{T}, Ce::Array{T}, 
                                r::T) where T

Intersection between the euclidan disks (Ce, r) and the euclidean line (P1P2)
"""
function intersect_linecircle(P1::Array{T}, P2::Array{T}, Ce::Array{T}, 
                              r::T) where T
    ε = T(10) ^(- 6)
    A, B, C =  P2[2] - P1[2], P1[1] - P2[1], P2[1] * P1[2] - P1[1] * P2[2]
    P, Q = if abs(B) > ε
        a = A ^ 2 + B ^ 2
        b = 2 * A * C + 2 * A * B * Ce[2] - 2 * B ^ 2 * Ce[1]
        c = C ^ 2 + 2 * B * C * Ce[2] - B ^ 2 * (r ^ 2 - Ce[1] ^ 2 - Ce[2] ^ 2)
        Δ = b ^ 2 - 4 * a * c 
        if Δ ≤ 0 
            @warn "there is no two intersection point calling  intersect_linecircle"
            return zeros(T, 2), zeros(T, 2)
        end
        x1, x2 = (- b - sqrt(Δ)) / (2 * a),  (- b + sqrt(Δ)) / (2 * a)  
        [x1, - (A * x1 + C) / B],  [x2, - (A * x2 + C) / B]
    else
        a = A ^ 2 + B ^ 2
        b = 2 * B * C + 2 * A * B * Ce[1] - 2 * A ^ 2 * Ce[2]
        c = C ^ 2 + 2 * A * C * Ce[1] - A ^ 2 * (r ^ 2 - Ce[1] ^ 2 - Ce[2] ^ 2)
        Δ = b ^ 2 - 4 * a * c 
        if Δ ≤ 0 
            @show Ce 
            @show r
            @show Δ
            @show A 
            @show B
            @show P1 
            @show P2
            @warn "there is no two intersection point calling  intersect_linecircle"
            return zeros(T, 2), zeros(T, 2)
        end
        y1, y2 = (- b - sqrt(Δ)) / (2 * a),  (- b + sqrt(Δ)) / (2 * a)  
        [- (B * y1 + C) / A, y1], [- (B * y2 + C) / A, y2]
    end
    return P, Q
end
#-------------------------------------------------------------------------------
"""
   
   I = intersect_lines(P1::Array{T}, P2::Array{T}, 
                       Q1::Array{T}, Q2::Array{T}) where T

Intersection of two lines given by points
"""
function intersect_lines(P1::Array{T}, P2::Array{T}, 
                         Q1::Array{T}, Q2::Array{T}) where T
    ε = 1e-8
    U, V = P2 - P1, Q2 - Q1
    if abs(U[1] * V[2] - U[2] * V[1]) < ε
        error("empty intersection in intersect_lines")
    end
    p = hcat(U, - V) \ (Q1 - P1)
    I = P1 + p[1] * U
    return I
end
#-------------------------------------------------------------------------------
function _map𝔻(p, c, rotθ::Complex = Complex(1))
    np = length(p)
    zp = [p[k][1] + im * p[k][2] for k = 1:np]
    iso_zp = [(z + c) / (1 + conj(c) * z) for z ∈ zp]
    iso_zp .*= rotθ
    return [[real(iso_zp[k]), imag(iso_zp[k])] for k = 1:np]
end
using Optim
function center𝔻(Plist::Vector{Vector{T}}) where T
    fcenter(c) = begin 
        zc = c[1] + c[2] * im 
        pc = _map𝔻(Plist, zc)
        #=
        minl = norm(pc[1] - pc[end])
        for k = 1:(length(pc) - 1)
            minl = min(minl, norm(pc[k] - pc[k + 1]))
            
        end
        return - minl 
        =#
        return maximum(norm.(pc))
    end
    res = optimize(fcenter, zeros(T, 2))
    @show res
    # output
    copt = Optim.minimizer(res)
    rotrandom = exp(2 * T(π) * rand(T) * im)
    return _map𝔻(Plist, copt[1] + copt[2] * im, rotrandom)
end
     
#-------------------------------------------------------------------------------
function plot_hexagon𝔻(Plist::Vector{Vector{T}}; 
                       debug::Bool = false, L::Any = nothing, 
                       prescribedl::Any = nothing) where T
    # recompute ploting data 
    if debug
        C4, r4, _ = perpandicular_by_P𝔻(Plist[3], Plist[4], Plist[4], T(0.1))
        b2 = C4[2]
        c2 = (r4 ^ 2 - norm(C4) ^ 2) / 2
        Corth = [0., (0.5 - c2) / b2]
        Rorth = sqrt(norm(Corth - C4) ^ 2 - r4 ^ 2)
        I1, I2 = Plist[(end - 1):end]
    end
    limxy, np = 1.5, 100
    if L == nothing
        _, L = Main.figure(1)
        L = L[1]
    end
    θ = LinRange(0, 2 * T(π), 200)
    # Poincare disk
    # colors
    cr = (rand(3)...,)
    cr = (0.2944, 0.3474, 0.6629)
    gr = (0.2, 0.8, 0.1)
    gb = (0.2, 0.1, 0.8)
    cred = (0.4, 0., 0.)
    cblue = (0., 0., 0.5)
    cgo = (0.3, 0.4, 0.2)
    sc = 5
    # plot a unite 𝔻 circle 
    Main.poly!(L, Circle(Point2f(0, 0.), 1.), color = RGBA(1., 0.1, 0.3, 0.3))
    Main.plot_points(hcat(cos.(θ), sin.(θ)), color = (0.8, 0.1, 0.1), 
                           scalep = 5, scene = L)
    # plot vertices
    Main.plot_points(Plist, color = cred, scalep = 15, scene = L)
    if debug
        # last geodesic 
        Main.poly!(L, Circle(Point2f(C4,), r4), color = RGBA(gr..., 0.2))
        Main.plot_points(vcat(I1', I2'), color = cr, scalep = 15, scene = L)
        # orth geodesic 
        Main.plot_points(Corth', color = cgo, scalep = 15, scene = L)
        Main.poly!(L, Circle(Point2f(Corth,), Rorth), color = RGBA(cgo..., 0.2))
    end
    # plot geodesics and  length
    lkall = zeros(6)
    # periodic indexing
    Ip(k) = mod(k - 1, length(Plist)) + 1
    for k = 1:length(Plist)
        #= Pk, Pkp = if k == length(Plist)
            Plist[1], Plist[end]
        else
           Plist[k], Plist[k + 1]
        end =#
        Pk, Pkp = Plist[Ip(k)],  Plist[Ip(k + 1)]
        C, rC, _, _, lk = geodesic𝔻(Pk, Pkp)
        if prescribedl != nothing
            @show lk
        end
        lkall[k] = lk
        pUV = arcparam𝔻(Pk, Pkp, C, rC, np)
        col = mod(k, 2) == 0 ? cred : cblue
        Main.plot_points(pUV, color = col, scalep = sc, scene = L)
    end 
    errororth = Inf
    for k = 1:length(Plist)
        Pk, Pkp = Plist[Ip(k)],  Plist[Ip(k + 1)]
        Pkm, Pk = Plist[Ip(k - 1)],  Plist[Ip(k)]
        Ckm, rkm, _ = geodesic𝔻(Pkm, Pk)
        Ckp, rkp, _ = geodesic𝔻(Pk, Pkp)
        # check orthogonality of consecutive geodesics 
        vorth =  abs(rkm ^ 2 + rkp ^ 2 - norm(Ckm - Ckp) ^ 2)
        if errororth > vorth 
            errororth = vorth
        end
    end
    # check orthogonality
    if prescribedl != nothing 
        error = maximum(abs.(lkall[1:2:end] - prescribedl))
        @printf "error               = %4.20f\n" error
        @printf "orthogonality error = %4.20f" errororth
    end
    Main.Makie.xlims!(L, -limxy, limxy)
    Main.Makie.ylims!(L, -limxy, limxy)
    return
end