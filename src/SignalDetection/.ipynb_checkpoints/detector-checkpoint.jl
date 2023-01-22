using Base: @kwdef
using ..DigitalModulation: Slicer, ModType


#=============== Matched Filter(MF) ===============#
struct MF end
function mf!(x::AbstractVector, y::AbstractVector, H::AbstractMatrix)
    x .= H' * y
    for i in eachindex(x)
        h_i = @view H[:,i]
        x[i] /= h_i'*h_i
    end
    return x
end


#=============== Zero Forcing(ZF) ===============#
struct ZF end
"""
ZF Detetction
"""
function zf(y::AbstractVector, H::AbstractMatrix)
    if size(H,1)==size(H,2)
        H \ y
    else
        (H'*H) \ (H'*y)
    end
end

function zf!(x::AbstractVector, y::AbstractVector, H::AbstractMatrix)
    if size(H,1)==size(H,2)
        x .= H \ y
    else
        x .= ((H'*H)^-1 * H') * y
    end
    x
end

"""
get ZF Weght Matrix
"""
function zf_weight(H::AbstractMatrix)
    if size(H,1)==size(H,2)
        return inv(H)'
    else
        return H*(H'H)^-1
    end
end

function zf!(x::AbstractVector, nvar::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0::Number)
    W = size(H,1)==size(H,2) ? H^-1 : (H'*H)^-1*H'
    nvar .= N0 .* slicedot(W', W')
    x .= W * y
    x, nvar
end

function zf!(x, nvar, y::AbstractVector, H::AbstractMatrix, N0::AbstractVector)
    W = conj(zf_weight(H))
    for i in eachindex(x)
        w = @views W[i,:]
        x[i] = dot(w, y)
        var = 0.0
        for j in eachindex(w)
            var += abs2(w[j]) * N0[j]
        end
        nvar[i] = var
    end
    x, nvar
end


function zf(y::AbstractVector{T}, H::AbstractMatrix{T}, N0::AbstractVector) where T
    x    = Array{T}(undef, size(H,2))
    nvar = Array{T}(undef, size(H,2))
    x, nvar = zf_soft!(x, nvar, y, H, N0)
    x, nvar
end


#=============== MMSE ===============#

struct MMSE end

"""
MMSE Detection
"""
mmse(y::AbstractVector, H::AbstractMatrix, N0::Number) = (H'*H + N0*I) \ (H'*y)
mmse(y::AbstractVector, H::AbstractMatrix, N0::AbstractVector) = H' * (H*H' + Diagonal(N0)) \ y

function mmse(y::AbstractVector, H::AbstractMatrix, N0::Number, α::Number)
    W = ((1+α^2)*(H*H') - 2α*H*H' + (size(y,1)*α + N0)*I)^-1 * H
    return W' * y
end

function mmse!(x̂, y, H, N0)
    x̂ .= mmse(y, H, N0)
    x̂
end
function mmse!(x̂, y, H, N0, α)
    x̂ .= mmse(y, H, N0, α)
    x̂
end

function mmse!(x̂, nvar, y, H, N0, x_prior)
    x̄, Σ = x_prior.mean, Diagonal(x_prior.var)
    nI = length(N0)==1 ? N0*I : Diagonal(N0)
    Wmmse = (H'*H + N0*I)^-1 * H'
    mu = real(slicedot(Wmmse', H))
    x̂ .= (Σ^-1 + H'*nI*H)^-1*H'*nI*y .+ Σ^-1 * x̄ # 事後推定
end


"""
MMSE Weight Matrix
"""
mmse_weight(H::AbstractMatrix, N0::Number) = H*inv(H'H + N0*I)
mmse_weight(H::AbstractMatrix, N0::AbstractVector) = inv(H*H' + Diagonal(N0))*H


# MMSE-soft
function mmse!(x̂::AbstractVector, nvar::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0::Number)
    W = H * (H'*H + N0*I)^-1
    μ = real(slicedot(W, H))
    x̂    .= (W'*y) ./ μ
    nvar .= real((1 .- μ) ./ μ)
    x̂, nvar
end
function mmse!(x̂::AbstractVector, nvar::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0::AbstractVector)
    W = (H*H' + Diagonal(N0))^-1 * H
    μ = real(slicedot(W, H))
    x̂ .= (W' * y) ./ μ
    nvar .= real((1 .- μ) ./ μ)
    return x̂, nvar
end

function mmse!(x̂::AbstractVector, nvar::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0::Number, α::Number)
    M = size(y,1)
    W = (1-α)*((1-α)^2*(H*H') + (M*α + N0)*I)^-1 * H
    μ = (1-α)*real(slicedot(W, H))
    x̂    .= (W'*y) ./ μ
    nvar .= real((1.0 .- μ) ./ μ)
end

function mmse!(x̂::AbstractVector, nvar::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0::AbstractVector, α::Number)
    M = size(y,1)
    W = (1-α)*((1-α)^2*(H*H') + (M*α .+ Diagonal(N0)))^-1 * H
    μ = (1-α)*real(slicedot(W, H))
    x̂    .= (W'*y) ./ μ
    nvar .= real((1.0 .- μ) ./ μ)
end

mmse(y::AbstractVector{T}, H::AbstractMatrix{T}, nI::AbstractMatrix) where T = return (H'*H + nI)\(H'*y)

function soft_mmse(y::AbstractVector{T}, H::AbstractMatrix{T}, nI::AbstractMatrix) where T
    x = Array{T}(undef, size(H,2))
    nvar = Array{Float64}(undef, size(H,2))
    soft_mmse!(x, y, H, nI)
    x, nvar
end

"""
    soft_mmse!(x, nvar, y, H, N0, xhat, vars)

SISO MMSE-PIC Detection
"""
function soft_mmse!(x, nvar, y::AbstractVector, H::AbstractMatrix, N0, xhat::AbstractVector, var::AbstractVector)
    Nt = length(x)
    r  = y .- H * xhat # 残差ベクトル
    nI = length(N0)==1 ? N0*I : Diagonal(N0)
    Δ  = Diagonal(var)
    Rinv = inv(H*Δ*H' + nI) # 残差共分散行列の逆行列

    for n = 1:Nt
        h = @view H[:,n]
        η = real(h' * Rinv * h)
        δ = 1 - Δ[n,n]
        γ = 1 / (1 + η*δ)
        μ = γ * η # 利得
        # @show μ
        # α = η / (1 + η*δ)
        # x[n] = (w'*r + α * xhat[n]) / μ
        x[n] = γ * (h'*Rinv*r + η*xhat[n]) / μ # 推定値
        ν = (1 - μ) / μ # 分散
        nvar[n] = ν
    end
    x, nvar
end



# 行列の各列ベクトル同士の内積をとる
function slicedot(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
   arr = Array{T}(undef, size(A,2))
   for j in axes(A,2)
       arr[j] = @views dot(A[:,j], B[:,j])
   end
   arr
end

# MMSE検出後の等価雑音分散を計算
function calc_equivar(w, y)
    μ = real(dot(w,y)) # 等価利得
    (1 - μ) / μ # 等価雑音分散
end



struct MRC end
"""
MRC Detection
"""
mrc(y, H) = H'y
mrc!(x, y, H) = (x .= mrc(y,H))

@kwdef struct SIC{T}
    detector::T = ZF()
    ordering::Symbol = :default
end
"""
Ordered SIC
"""
function sic!(x, y::AbstractArray{T}, H::AbstractArray{T}, N0, slicer, sort=:SIR, criterion=:SIR) where T
    nrx, ntx = size(H)
    r = copy(y)
    nI = length(N0)==1 ? N0*I : Diagonal(N0) # 雑音共分散行列
    order = columnsort(H, sort, rev=true)
    @views _sic!(x, r, H, nI, order, slicer) #SICの実行
    return x
end

# sicの内部関数
function _sic!(x, r, H, nI, order, slicer)
    R = H*H' + nI
    w = similar(r)
    for (st,i) in enumerate(order)
        hi = @view H[:,i]
        mul!(w, inv(R), hi)
        tmp = dot(w, r)
        x[i] = slicer(tmp)
        # symbol subtraction
        r .-= hi * x[i]
        st < length(x) && (R .-= hi*hi')
    end
end
function qrsic!(x, y, H, slicer)
    nrx, ntx = size(H)
    ord = columnsort(H, :SIR, rev=false)
    Q, R = @views qr(H[:,ord])
    r = (Q'y)[1:ntx]
    @views _qrsic!(x[ord], r, R, slicer)
end

function _qrsic!(x, r, R, slicer)
    ntx = size(R,2)
    for n = ntx:-1:1
        tmp = r[n]
        for m = n+1:ntx # SIC
            tmp -= R[n,m] * x[m]
        end
        tmp /= R[n,n]
        x[n] = slicer(tmp)
    end
end

struct VBLAST{modtype}
    detection::Symbol
    slicer::Slicer
end
function VBLAST(m::ModType, detection=:ZF)
    slicer = Slicer(m)
    VBLAST{typeof(m)}(detection, slicer)
end
"""
    vblast!(x, y, H, slicer) -> ZF
    vblast!(x, y, H, N0, slicer) -> MMSE

V-BLAST Detection
    ZF : argmin_{i} |w_i|
    MMSE: argmin_{i} mu_i - mu_i^2
"""
# MMASEによる検出
function vblast!(x::AbstractVector, y::AbstractVector, H::AbstractMatrix, N0, slicer, detection)
    if detection == :ZF
        x = @views _vblast_by_zf!(x, copy(y), H, slicer) #SICの実行
    else
        x = @views _vblast_by_mmse!(x, copy(y), H, N0, slicer) #の実行
    end
    return x
end

vblast!(x, y, H, N0; slicer, detection) = vblast!(x, y, H, N0, slicer, detection)


# V-BLAST by ZF weight
function _vblast_by_zf!(x, r, H, slicer)
    idxs = [1:length(x);]
    for i = eachindex(x)
        W = zf_weight(H[:,idxs])
        if i < length(x)
            k = map(norm, eachcol(W)) |> argmin
        else
            k = 1
        end
        idx = idxs[k] # 最小ノルムを選択
        x[idx] = slicer(dot(W[:,k], r)) # 判定
        r .-= H[:,idx] * x[idx]
        deleteat!(idxs, k) # 選択したindexを候補から除外
    end
    x
end


# V-BLAST by MMSE weight
function _vblast_by_mmse!(x, r, H, N0, slicer)
    idxs = [1:length(x);]
    for i = 1:length(x)
        W = mmse_weight(H[:,idxs], N0)
        k = map(norm, eachcol(W)) |> argmin # SNR最大化フィルタを選択
        # k = @views _maximum_post_sinr(W, H[:,idxs])
        idx = idxs[k] # 最小ノルムを選択
        x[idx] = slicer(dot(W[:,k], r)) # 判定
        r .-= H[:,idx] * x[idx]
        deleteat!(idxs,k) # 選択したindexを候補から除外
    end
    x
end

function _maximum_post_sinr(W, H)
    mu = slicedot(W,H)
    var = real(mu .- abs2.(mu))
    return argmin(var)
end



"""
PIC Detection
"""
function pic!(x, nvar, y, H, N0, mod_type, n_iter; LLRs=nothing, σs2=1.0) where T
    m, n = size(H)

    # 初期推定
    if isnothing(LLRs)
        soft_mmse!(x, nvar, y, H, N0)
        LLRs = demod(mod_type, x, nvar)
    end

    # PICの実行
    r    = similar(y)              # 干渉除去信号
    wh   = zeros(T, 1, n)          # ウェイトベクトル
    x_soft    = zeros(T, n)        # ソフトシンボルベクトル
    Rinv = zeros(T, m, n)          # ウェイトベクトル計算のための逆行列
    σe2  = zeros(Float64, n)       # 推定分散
    δ    = zeros(Float64, n)       # 推定誤差分散ベクトル
    nI = length(N0)==1 ? N0*I : Diagonal(N0) # 雑音分散行列
    iter = 0
    while iter < n_iter
        x_soft = mod(mod_type, LLRs)
        s .= getmean.(soft_data)
        δ .= getvar.(soft_data)
        _pic!(x, nvar, y, H, nI, s, r, wh, Rinv, σs2, δ)   # PICの実行
        iter += 1
    end
end

# 並列干渉除去
function _pic!(x, nvar, y, H, nI, x_soft, r, wh, Rinv, σ2s, δ)
    Δ = Diagonal(δ)             # 対角行列
    Rinv .= inv(H*Δ*H' + nI)    # 残留ノイズの精度行列

    # PICの実行
    γ = μ = 0.0
    for i in axes(H,2)
        # 干渉除去
        r .= @views y .- H[:,Not(i)] * x_soft[Not(i)]

        # MMSEウェイトベクトルwiの計算
        h = @view H[:,i]                      # チャネルベクトル
        γ  = real(h'*Rinv*h)                  # γ
        μ  = γ / (1.0 + γ*(σs2-δ[i]))         # 利得
        wh .= h'*Rinv * (1.0 / (1.0 + γ*(σ2s-δ[i])))  # 重み行ベクトル

        # 推定量と分散の計算
        x[i] = (wh * r / μ)[1]  # 推定値
        nvar[i]  = (1 - μ) / μ  # 雑音分散
    end
    x, nvar
end
