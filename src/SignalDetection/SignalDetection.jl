# SigDet.jl
module SignalDetection

# source files
include("./detector.jl")
include("./mld.jl")     # MIMO-MLD 
# include("./qrm_mld.jl") # 
include("./MIMO/MLD/SphereDetector.jl")



# packages
import Base: @kwdef
using LinearAlgebra: Diagonal, diagm, I
using ..DigitalModulation
using .SphereDetector: sphere_detector

to_real(x::AbstractVector{T}) where T<:Complex = [real(x); imag(x)]
to_real(A::AbstractMatrix{T}) where T<:Complex = [real(A) -imag(A); imag(A) real(A)]
function to_complex(x::AbstractVector{T}) where T<:Real
    n = length(x)
    return x[1:div(n,2)] + im*x[div(n,2)+1:end]
end

export signal_detection, detect!, detect
export MF, ZF, MRC, MMSE, SIC, VBLAST, PIC, MLD
export to_real, to_complex

# MF detector
detect!(::MF, x::AbstractVector, y::AbstractVector, H::AbstractMatrix, args...) = mf!(x, y, H)


# ZF detector
detect!(::ZF, x::AbstractVector, y::AbstractVector, H::AbstractMatrix, args...) = zf!(x, y, H)
detect!(::ZF, x::AbstractVector, nvars::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0, args...) = zf!(x, nvars, y, H, N0)

# detect!(x::AbstractVector, y::AbstractVector, H::AbstractMatrix, args...; method::ZF) = zf!(x, y, H)
(::ZF)(y::Number,h::Number; kwargs...) = zf(y,h)
(::ZF)(x::T, y::T,h::T; kwargs...) where T <: AbstractVector = zf(y,h)
(::ZF)(x::AbstractVector,y::AbstractVector,H::AbstractMatrix; kwargs...) = zf!(x, y, H)


# MRC detector
detect!(::MRC, x::AbstractVector, y::AbstractVector, h::AbstractVector, kwargs...) = mrc!(x, y, h)

# MMSE detector
detect(::MMSE, y::AbstractVector, H::AbstractMatrix, N0; kwargs...) = mmse(y, H, N0)
detect!(::MMSE, x::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0; kwargs...) = mmse!(x, y, H, N0)
detect!(::MMSE, x::AbstractVector, nvars::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0; kwargs...) = mmse!(x, nvars, y, H, N0)


(::MMSE)(y::T, h::T, N0; kwargs...) where T<:Number = mmse(y, h, N0)
(::MMSE)(y::AbstractVector{T}, H::AbstractMatrix{T}, N0; kwargs...) where T  = mmse(y, H, N0)
function (::MMSE)(x::AbstractVector{T}, y::AbstractVector{T}, H::AbstractMatrix{T}, N0; kwargs...) where T
    α = get(kwargs, :α, 0.0)
    if α == 0.0
        mmse!(x, y, H, N0)
    else
        mmse!(x, y, H, N0, α)
    end
end
function (::MMSE)(x::AbstractVector{T}, nvar, y::AbstractVector{T}, H::AbstractMatrix{T}, N0; kwargs...) where T
    α = get(kwargs, :α, 0.0)
    if α == 0.0
        x, nvar = soft_mmse!(x, nvar, y, H, N0)
    else
        x, nvar = soft_mmse!(x, nvar, y, H, N0, α)
    end
end
function (::MMSE)(x::AbstractVector, nvar::AbstractVector, y::AbstractVector{T}, H::AbstractMatrix{T}, N0, xhat; kwargs...) where T
    x, nvar = soft_mmse!(x, nvar, y, H, N0, getmean.(xhat), getvar.(xhat))
end

function detect!(::MMSE, x, y, H, nvar, args...)
    if args[:output] == :hard
        mmse!(x, y, H, nvar)
    else
        mmse!(x, y, H, nvar)
    end
end


(v::VBLAST)(x, y, H; kwargs...) = vblast!(x, y, H, kwargs[:N0], v.slicer, v.detection)
(v::VBLAST)(x, y, H, N0; kwargs...) = vblast!(x, y, H, N0, v.slicer, v.detection)


Base.show(io::IO, m::MIME"text/plain", mld::MLD) = print(io, "MIMO ML Detector")

detect(ml::MLD, y, H; kwargs...) = mld(y, H, ml.refs) 
(mld::MLD)(y, H; kwargs...) = mld(y, H, mld.refs)
(mld::MLD)(x, y, H; kwargs...) = mld!(x, y, H, mld.refs)



struct QRMLD{T}
    refs::Vector{T}
    K::Int
end
"""
    QRMLD(modtype; K)

"""
function QRMLD(modtype; K)
    refs = make_refs(modtype)
    QRMLD{eltype(refs)}(refs, K)
end
function (qrmld::QRMLD)(y, H; kwargs...)
    R = get(kwargs, :R, nothing)
    if isnothing(R)
        x, list = qrmld(y, H, qrmld.refs, qrmld.K)
    else
        order = get(kwargs, :order, nothing)
        x, list = qrmld(y, R, order, qrmld.refs, qrmld.K)
    end
    x
end
function (qrmld::QRMLD)(x, y, H; kwargs...)
    R = get(kwargs, :R, nothing)
    if isnothing(R)
        x, list = qrmld!(x, y, H, qrmld.refs, qrmld.K)
    else
        order = get(kwargs, :order, nothing)
        x, list = qrmld!(x, y, R, order, qrmld.refs, qrmld.K)
    end
    x
end

struct PIC end

struct Dim{N} end
Dim(N) = Dim{N}()

is_nvar_required(detector) = detector isa Union{MMSE,SIC,VBLAST,PIC}

#=========================== Interface functions =============================#

"""
    signal_detection(y::AbstractArray, H::AbstractArray, nvar=nothing, method="ZF", option=())

信号検出方式
    SISO: :ZF, :MMSE
    SIMO: :MRC
    MIMO: :ZF, :MMSE,
"""
function signal_detection(y, H, N0=nothing; dims=1, detector=ZF(), output_type=:hard, xhat=nothing, kwargs...)
    # データの並び次元
    if dims==1
        nsym, nrx, ntx = size(y,1), size(H,2), size(H,3)
    elseif dims == 2
        nsym, nrx, ntx = size(y,2), size(y,1), size(H,1)
    end
    # 推定シンボル配列，等価雑音配列の確保
    x = Array{eltype(y)}(undef, ntx, nsym)
    if output_type == :hard && !is_nvar_required(detector)
        N0 = nothing
    end
    if nrx == ntx == 1
        x = vec(x); y = vec(y); H = vec(H)
    end
    if output_type==:hard
        if isnothing(xhat)
            routine!(detector, x, y, H, N0, Dim(dims); kwargs...)
        else
            routine!(detector, x, y, H, N0, Dim(dims), xhat; kwargs...)
        end
    else
        nvar = similar(x, Float64)
        if isnothing(xhat)
            routine!(detector, x, nvar, y, H, N0, Dim(dims); kwargs...)
        else
            routine!(detector, x, nvar, y, H, N0, Dim(dims), xhat; kwargs...)
        end
        return x, nvar
    end
end

function split_string(str)
    sym = split(uppercase(str), "_")
    if length(sym)==1
        return Symbol(sym[1]), nothing
    else
        return Symbol(sym[1]), parse(Int,sym[2])
    end
end

#SISO
function routine!(f, x::AbstractArray{T,1}, y::AbstractArray{T,1}, h::AbstractArray{T,1}, N0, ::Dim{1}; kwargs...) where T
    if N0 isa Nothing
        for i in eachindex(y)
            x[i] = f(y[i], h[i]; kwargs...)
        end
    elseif N0 isa Number
        for i in eachindex(y)
            x[i] = f(y[i], h[i], N0; kwargs...)
        end
    elseif N0 isa Vector
        for i in eachindex(y)
            x[i] = f(y[i], h[i], N0[i]; kwargs...)
        end
    else
        error()
    end
    x
end

# SIMO/MIMO
function routine!(m, x::AbstractArray{T,2}, y::AbstractArray{T,2}, H::AbstractArray{T,3}, N0, ::Dim{1}; kwargs...) where T
    if N0 isa Nothing
        for i in axes(y,1)
            @views detect!(m, x[:,i], y[i,:], H[i,:,:]; kwargs...)
        end
    elseif N0 isa Number
        for i in axes(y,1)
            @views detect!(m, x[:,i], y[i,:], H[i,:,:], N0; kwargs...)
        end
    elseif ndims(N0)==1
        for i in axes(y,1)
            @views detect!(m, x[:,i], y[i,:], H[i,:,:], N0; kwargs...)
        end
    elseif ndims(N0)==2
        for i in axes(y,1)
            @views detect!(m, x[:,i], y[i,:], H[i,:,:], N0[i,:]; kwargs...)
        end
    end
    x
end
function routine!(f, x::AbstractArray{T,2}, y::AbstractArray{T,2}, H::AbstractArray{T,3}, N0, ::Dim{2}; kwargs...) where T
    if N0 isa Nothing
        for i in axes(y,2)
            x[:,i] = @views detect(f, y[:,i], H[:,:,i]; kwargs...)
        end
    elseif N0 isa Number
        for i in axes(y,2)
            @views detect!(f, x[:,i], y[:,i], H[:,:,i], N0; kwargs...)
        end
    elseif ndims(N0)==1
        for i in axes(y,2)
            @views detect!(f, x[:,i], y[:,i], H[:,:,i], N0; kwargs...)
        end
    elseif ndims(N0)==2
        for i in axes(y,2)
            @views detect!(f, x[:,i], y[i,:], H[:,:,i], N0[i,:]; kwargs...)
        end
    end
    x
end

function routine!(f, x::AbstractArray{T,2}, nvar::AbstractArray{U,2}, y::AbstractArray{T,2}, H::AbstractArray{T,3}, N0, ::Dim{1}; kwargs...) where {T,U}
    if N0 isa Nothing
        for i in axes(y,1)
            @views detect!(f, x[:,i], nvar[:,i], y[i,:], H[i,:,:]; kwargs...)
        end
    elseif ndims(N0)<=1
        for i in axes(y,1)
            @views detect!(f, x[:,i], nvar[:,i], y[i,:], H[i,:,:], N0; kwargs...)
        end
    elseif ndims(N0)==2
        for i in axes(y,1)
            @views detect!(f, x[:,i], nvar[:,i], y[i,:], H[i,:,:], N0[i,:]; kwargs...)
        end
    end
    x, nvar
end

function routine!(f, x::AbstractArray{T,2}, nvar::AbstractArray{U,2}, y::AbstractArray{T,2}, H::AbstractArray{T,3}, N0, ::Dim{2}; kwargs...) where {T,U}
    if N0 isa Nothing
        for i in axes(y,dims)
            @views f(x[:,i], nvar[:,i], y[:,i], H[:,:,i]; kwargs...)
        end
    elseif N0 isa Number
        for i in axes(y,2)
            @views f(x[:,i], nvar[:,i], y[:,i], H[:,:,i], N0; kwargs...)
        end
    elseif ndims(N0)==1
        for i in axes(y,2)
            @views f(x[:,i], nvar[:,i], y[:,i], H[:,:,i], N0; kwargs...)
        end
    elseif ndims(N0)==2
        for i in axes(y,2)
            @views f(x[:,i], nvar[:,i], y[:,i], H[:,:,i], N0[i,:]; kwargs...)
        end
    end
    x, nvar
end
function routine!(f, x::AbstractArray{T,2}, nvar::AbstractArray{U,2}, y::AbstractArray{T,2}, H::AbstractArray{T,3}, N0, ::Dim{1}, xhat::AbstractArray{V,2}; kwargs...) where {T,U,V}
    if N0 isa Nothing
        for i in axes(y,dims)
            @views f(x[:,i], nvar[:,i], y[i,:], H[i,:,:]; kwargs...)
        end
    elseif N0 isa Number
        for i in axes(y,1)
            @views f(x[:,i], nvar[:,i], y[i,:], H[i,:,:], N0, xhat[i,:]; kwargs...)
        end
    elseif ndims(N0)==1
        for i in axes(y,1)
            @views f(x[:,i], nvar[:,i], y[i,:], H[i,:,:], N0, xhat[i,:]; kwargs...)
        end
    elseif ndims(N0)==2
        for i in axes(y,1)
            @views f(x[:,i], nvar[:,i], y[i,:], H[i,:,:], N0[i,:]; kwargs...)
        end
    end
    x, nvar
end
function routine!(f, x::AbstractArray{T,2}, nvar::AbstractArray{U,2}, y::AbstractArray{T,2}, H::AbstractArray{T,3}, N0, ::Dim{2}, xhat::AbstractArray{V,2}; kwargs...) where {T,U,V}
    if N0 isa Nothing
        for i in axes(y,dims)
            @views f(x[:,i], nvar[:,i], y[:,i], H[i,:,:]; kwargs...)
        end
    elseif N0 isa Number
        for i in axes(y,2)
            @views f(x[:,i], nvar[:,i], y[:,i], H[:,:,i], N0, xhat[i,:]; kwargs...)
        end
    elseif ndims(N0)==1
        for i in axes(y,2)
            @views f(x[:,i], nvar[:,i], y[:,i], H[:,:,i], N0, xhat[i,:]; kwargs...)
        end
    elseif ndims(N0)==2
        for i in axes(y,2)
            @views f(x[:,i], nvar[:,i], y[:,i], H[:,:,i], N0[i,:]; kwargs...)
        end
    end
    x, nvar
end

#=

"""
    siso_zf(y, h, nvar, output)

SISO-ZF Detection
"""
function siso_zf(y, h)
    y = vec(y)
    h = vec(h)
    x = similar(y)
    f = zf
    routine!(f, x, y, h)
    return x
end

function siso_zf_soft(y, h, N0)
    y  = vec(y)
    h  = vec(h)
    N0 = length(N0)==1 ? N0 : vec(N0)
    x    = similar(y)
    nvar = similar(y, Float64)
    routine!(zf_soft, x, nvar, y, h, N0)
    return x, nvar
end

"""
    siso_mmse(y, h, nvar, output)

SISO-MMSE Detection
"""
function siso_mmse(y, h, nvar, output)
    if output=="hard"
        return mmse(y, h, nvar)
    else
        return mmse_soft(y, h, nvar)
    end
end


"""
    simo_mrc(y, h, nvar, output)

MRC(Maximum-Ratio-Combining)
"""
function simo_mrc(y, h, N0)
    x̂ = Array{eltype(y)}(undef, size(y,2)); # 出力配列
    for n in axes(y,2)
        x̂[n] = @views mrc(y[:,n], h[:,1,n])
    end
    x̂
end


"""
MIMO ZF-Detection
"""
function mimo_zf(y, H)
    n_rx, n_tx, n_data = size(H,1), size(H,2), size(y,2)
    x = Array{eltype(y)}(undef, n_tx, n_data)

    if ndims(H) == 2
        for n in axes(y,2)
            x[:,n] = @views zf(y[:,n], H)
        end
    else
        for n in axes(y,2)
            x[:,n] = @views zf(y[:,n], H[:,:,n])
        end
    end
    x
end

function mimo_zf_soft(y, H, N0::Number)
    n_rx, n_tx, n_data_set = size(H,1), size(H,2), size(y,2)
    x_hat = Array{eltype(y)}(undef, n_tx, n_data_set)
    eq_nvar = Array{Float64}(undef, n_tx, n_data_set)

    for n in axes(y,2)
        @views zf_soft!(x_hat[:,n], eq_nvar[:,n], y[:,n], H[:,:,n], N0)
    end
    x_hat, eq_nvar
end

function mimo_zf_soft(y, H, N0::AbstractVector)
    n_rx, n_tx, n_data_set = size(H,1), size(H,2), size(y,2)
    x_hat = Array{eltype(y)}(undef, n_tx, n_data_set)
    eq_nvar = Array{Float64}(undef, n_tx, n_data_set)

    for n in axes(y,2)
        @views zf_soft!(x_hat[:,n], eq_nvar[:,n], y[:,n], H[:,:,n], N0)
    end
    x_hat, eq_nvar
end

function mimo_zf_soft(y, H, N0::AbstractMatrix)
    n_rx, n_tx, n_data_set = size(H,1), size(H,2), size(y,2)
    x_hat = Array{eltype(y)}(undef, n_tx, n_data_set)
    eq_nvar = Array{Float64}(undef, n_tx, n_data_set)

    for n in axes(y,2)
        @views zf_soft!(x_hat[:,n], eq_nvar[:,n], y[:,n], H[:,:,n], N0[:,n])
    end
    x_hat, eq_nvar
end


"""
    mimo_mmse(y, H, N0::Number)
    mimo_mmse(y, H, N0::AbstractVector)
    mimo_mmse(y, H, N0::AbstractMatrix)

MIMO-MMSE-detection
"""
function mimo_mmse(y, H, N0::Number)
    x = Array{eltype(y)}(undef, size(H,2), size(y,2)); # 信号推定配列
    for n in axes(y,2)
        @views mmse!(x[:,n], y[:,n], H[:,:,n], N0*I)
    end
    x
end

function mimo_mmse(y, H, N0::AbstractVector)
    x = Array{eltype(y)}(undef, size(H,2), size(y,2))
    for n in axes(y,2)
        @views mmse!(x[:,n], y[:,n], H[:,:,n], N0)
    end
    x
end

function mimo_mmse(y, H, N0::AbstractMatrix)
    x = Array{eltype(y)}(undef, size(H,2), size(y,2));
    for n in axes(y,2)
        @views mmse!(x[:,n], y[:,n], H[:,:,n], N0[:,n])
    end
    x
end

"""
    mimo_mmse_soft(y, H, N0::Number)

"""
function mimo_mmse_soft(y, H, N0, ntx, nrx, dims)
    x = Array{eltype(y)}(undef, ntx, size(y,dims))
    nvar  = Array{Float64}(undef, ntx, size(y,dims))

    routine!(mmse_soft!, x, nvar, y, H, N0, dims)
    x, nvar
end
#=
function mimo_mmse_soft(y, H, N0::AbstractVector)
    x = Array{eltype(y)}(undef, size(H,2), size(y,2))
    nvar  = Array{Float64}(undef, size(H,2), size(y,2))
    for n in axes(y,2)
        @views mmse_soft!(x[:,n], nvar[:,n], y[:,n], H[:,:,n], N0)
    end
    x, nvar
end
function mimo_mmse_soft(y, H, N0::AbstractMatrix)
    x     = Array{eltype(y)}(undef, size(H,2), size(y,2));
    nvar  = Array{Float64}(undef, size(H,2), size(y,2))
    for n in axes(y,2)
        @views mmse_soft!(x[:,n], nvar[:,n], y[:,n], H[:,:,n], N0[:,n])
    end
    x, nvar
end
=#

"""
    mimo_sic(y, H, N0, mod_type)
MIMO SIC-Detection
"""
function mimo_sic(y, H, N0, args) where T
    slicer = x -> slice(mod_type, x)
    routine!(sic!, x, y, H, N0; kwargs=args)
    x
end


function mimo_sic_soft(y, H::AbstractMatrix{T}, N0::Number, mod_type) where T
    x_hat = Array{eltype(y)}(undef, size(H,2), size(y,2)) # 推定ベクトル列
    nvar  = similar(x_hat)                                # 等価雑音分散

    for n in axes(y,2)
        @views sic_soft(x_hat[:,n], nvar[:,n], y[:,n], H,  N0*I)
    end
    x_hat, nvar
end

function mimo_sic_soft(y, H::AbstractArray{T,3}, N0::Number, mod_type) where T
    x_hat = Array{eltype(y)}(undef, size(H,2), size(y,2)) # 推定ベクトル列
    nvar  = similar(x_hat)                                # 等価雑音分散

    for n in axes(y,2)
        @views sic_soft(x_hat[:,n], nvar[:,n], y[:,n], H[:,:,n],  N0*I)
    end
    x_hat, nvar
end

"""
    mimo_vblast(y, H, N0, mod_type, weight=:MMSE)
VBLAST検出
"""
function mimo_vblast(y, H, N0, mod_type)
    x = Array{eltype(y)}(undef, size(H,2), size(y,2)) # 推定ベクトル列
    for n in axes(y,2)
        @views vblast!(x[:,n], y[:,n], H[:,:,n], N0, mod_type, :MMSE)
    end
    x
end



"""
    mimo_pic(y,H,N0,output,qam,iter)

MIMO-PIC(Parallel Interference Cancellation)
"""
function mimo_pic(y, H, N0::Number; mod_type, n_iter)
    ntx, ntx = size(H)[1:2]
    x̂ = Array{eltype(y)}(undef, ntx, size(y,2))
    eq_nvar = Array{Float64}(undef, size(x̂))
    for n in axes(y,2)
        @views pic!(x̂[:,n], eq_nvar[:,n], y[:,n], H[:,:,n], N0, mod_type, n_iter)
    end
    return x̂, eq_nvar
end

function mimo_pic(y, H, N0::AbstractVector; mod_type, n_iter)
    Nr, Nt = size(H)[1:2]
    x = Array{eltype(y)}(undef, Nt, size(y,2))
    eq_nvar = Array{Float64}(undef, size(x_hat))
    for n in axes(y,2)
        @views pic!(xhat[:,n], eq_nvar[:,n], y[:,n],H[:,:,n], N0, mod_type, n_iter)
    end
    return x_hat, eq_nvar
end

function mimo_pic(y, H, N0::AbstractMatrix; mod_type, n_iter)
    Nr, Nt = size(H)[1:2]
    x = Array{eltype(y)}(undef, Nt, size(y,2))
    eq_nvar = Array{Float64}(undef, size(x_hat))
    for n in axes(y,2)
        @views pic!(xhat[:,n], eq_nvar[:,n], y[:,n], H[:,:,n], N0[:,n], mod_type, n_iter)
    end
    return x_hat, eq_nvar
end


"""
    mimo_mmse_pic(y, H, N0, xp, evar, output)

MIMO-MMSE-PIC(Succesive Interference Cancellation) detection
"""
function mimo_mmse_pic(y, H, N0, xp, evar, output)
    """
    mimo_mmse_pic: MMSE-PIC検出
    """
    Nr, Nt = size(H) # 送受信アンテナ数
    x_hat = Array{eltype(y)}(undef, Nt, size(y,2)) # 推定信号配列
    nI = Array{eltype(y)}(N0[1]*I,Nt,Nt) # 雑音分散行列
    if output=="hard"
        for n in 1:size(y,2)
            if length(N0)>1
                nI[diagind(nI)] .= view(N0,:,n)
            end
            x_hat[:,n] = mmse_pic(view(y,:,n),view(H,:,:,n),nI,view(xp,:,n))[1]
        end
        return x_hat
    elseif output==:soft
        eq_nvar = Array{Float64}(undef, size(x_hat))
        for n in 1:size(y,2)
            if length(N0)>1
                nI[diagind(nI)] .= view(N0,:,n)
            end
            x_hat[:,n], eq_nvar[:,n] = pic_detection(view(y,:,n),view(H,:,:,n),nI,view(xp,:,n),view(evar,:,n))
        end
        return x_hat, eq_nvar
    end
end


"""
    mimo_mld()

MIMO-MLD: Maximum Likelihood Detection
"""
function mimo_mld(y, H, mod_type)
    ntx, ntx, ndata = size(H)
    x = Array{eltype(y)}(undef, ntx, ndata)
    refs = getrefs(mod_type)
    for n in axes(y,2)
        @views mld!(x[:,n], y[:,n], H[:,:,n], refs)
    end
    return x
end

"""

MIMO-QRSIC: QR分解によるSIC検出
"""
function mimo_qrsic(y, H, mod_type)
    ntx, ntx, ndata = size(H)
    x = Array{eltype(y)}(undef, ntx, ndata)
    slicer = x -> slice(mod_type, x)
    for n in axes(y,2)
        @views qrsic!(x[:,n], y[:,n], H[:,:,n], slicer)
    end
    return x
end

"""

QR-MLD
    argments:
        refs : referene symbols
        M    : number of survive paths for M algorithm

"""
function mimo_qrmld(y::AbstractArray{T}, H::AbstractArray{T}, mod_type, M) where T
    x = Array{T}(undef, size(H,2), size(H,3)); # 出力配列
    refs = getrefs(mod_type)
    for n in axes(y,2)
        x[:,n] .= @views qrmld(y[:,n], H[:,:,n], refs, M)[1]
    end
    return x
end

"""
    mimo_sdmld()

MIMO-SDMLD: Sphere Decoding ML Detector
"""
function mimo_sdmld(y::AbstractArray{T}, H::AbstractArray{T}, mod_type, K) where T
    x = Array{T}(undef, size(H,2), size(H,3)); # 出力配列

    for n in axes(y,2)
        @views sdmld!(x[:,n], y[:,n], H[:,:,n], mod_type, K)
    end
    x
end



=#


end # module
