# cir_est.jl : インパルス応答(CIR)の推定


"""
実際の連続応答からシンボルブロックの平均CIR応答を算出する
"""
function resize(cir::AbstractArray{<:Complex,4}, fft_size, cp_size=0; average=false, kwargs...)
    L = size(cir,1) # チャネル長

    # シンボルブロックに分割
    cir = reshape(cir, size(cir,1), fft_size+cp_size, :, size(cir,3), size(cir,4))
    cir = cir[:, cp_size+1:end,:,:,:] # cp部分除去
    cir = dropdims(sum(cir, dims=2) ./ fft_size, dims=2); # サンプル時間平均

    # 時間方向平均
    if average
        f_size = size(cir, 2)
        w_size = get(kwargs, :w_size, f_size)
        w_size = ifelse(w_size>f_size, f_size, w_size)
        divs = cld(f_size, w_size)
        if mod(f_size, w_size) > 0 # 割り切れない場合
            ave_cir = @views [
                mean(cir[:,i:i+w_size-1,:,:], dims=2) for i = 1:w_size:w_size*(divs-1)
            ]
            @views push!(ave_cir, mean(cir[:,w_size*(divs-1):end, dims=2]))
        else
            ave_cir = @views [
                mean(cir[:,i:i+w_size-1,:,:], dims=2) for i = 1:w_size:f_size
            ]
        end
        ave_cir = cat(ave_cir..., dims=2)
        cir = ave_cir
    end
    return cir
end
resize(ofdm::Ofdm , cir::AbstractArray{<:Complex,4}; average=false, divs=1) =
    resize(cir, ofdm.n_fft, ofdm.n_gi; average, divs)

"""
CFRをCIRに変換する
"""
to_cir(ofdm::Ofdm, CFR, L=size(CIR,1)-1) = return fft(CFR, 1)[1:L+1, :, :]


"""
    cirest()
CIR(Channel Impulse Response)推定
    Inputs:
            L : max path length (<= Pilot Length)

    Outputs:
            x : estimated CIR
            idxs(option) : estimated nonzero path indices
"""
function cirest(ofdm::Ofdm, Y::AbstractArray{T}; L=ofdm.n_gi, pidxs=ofdm.idxs[:pilot], kwargs...) where T
    method = get(kwargs, :method, :LS)
    Φp = get_dft_matrix(ofdm, L)
    if method != :interp_LS
        Yp = @views Y[pidxs,:,:]
    else
        Yp = Y
    end
    cirest(Yp, Φp; kwargs...)
end
function cirest(Yp::AbstractArray, Φp::AbstractMatrix; nvar=nothing, idxs=Int[], method=:LS, args=(;), kwargs...)
    nrx, ntx = size(Yp,2), size(Yp,3)
    h = zeros(eltype(Yp), size(Φp,2), nrx, ntx) # CIR vectors

    if method == :LS
        h = _ls!(h, Yp, Φp, idxs)
    elseif method == :LMMSE
        h = _lmmse!(h, Yp, Φp, nvar; args...)
    elseif method == :OMP
        h, idxs = _omp!(h, Yp, Φp, nvar, idxs; args...)
    elseif method == :SOMP
        h, idxs = _somp!(h, Yp, Φp, nvar, idxs; args...)
    elseif method == :SBL
        h, idxs = _sbl!(h, Yp, Φp, nvar, idxs; args...)
    elseif method == :BMP
        h, idxs = _bmp!(h, Yp, Φp, nvar, idxs; args...)
    elseif method == :interp_LS
        h = _interp_ls!(h, Yp, nvar; args...);
    else
        error("methodが正しくありません。")
    end
    # @show idxs
    return reshape(h, :, nrx, ntx), idxs
end

# ls!のラッパー
function _ls!(x, y, Φ; idxs=nothing)
    x = reshape(x, size(x,1), :)
    y = reshape(y, size(y,1), :)
    idxs = isempty(idxs) ? nothing : idxs
    if isnothing(idxs)
        ls!(x, y, Φ)
    else
        ls!(x, y, Φ, idxs)
    end
    x
end

# lmmse!のラッパー
function _lmmse!(x, y, Φ, nvar; idxs=nothing, Σ=nothing, prior=nothing, kwargs...)
    if nvar isa Number
        x = reshape(x, size(x,1), :); y = reshape(y, size(y,1), :)
        x = lmmse!(x, y, Φ, nvar, idxs, Σ)
    else
        for j = axes(y,2)
            @views lmmse!(x[:,j,:], y[:,j,:], Φ, nvar[j], idxs, Σ)
        end
    end
    x
end


# OMPのラッパー関数
function _omp!(x, y, Φ, nvar, idxs; args...)
    n_rx = size(y,3)
    n_tx = size(x,3)

    if size(y,2) == size(y,3) == 1
        x = vec(x); y = vec(y)
    else
        x = reshape(x, size(x,1), :)
        y = reshape(y, size(y,1), :)
    end
    max_iter = get(args, :max_iter, size(y,1))
    max_iter = max_iter > size(y,1) ? size(y,1) : max_iter
    n = size(y,1)
    for i = 1:n_rx
        if ndims(nvar)==1
            ϵ = nvar * sqrt(n*l + 2*sqrt(n*log(n)))
        else
            ϵ = sqrt(nvar[i]) * sqrt(n*l + 2*sqrt(n*log(n)))
        end
        for j = 1:n_tx
            _, S = @views Greedy.omp!(x[:,2(i-1)+j], y[:,2(i-1)+j], Φ; ϵ=ϵ, max_iter=max_iter)
        end
        union!(idxs, S)
    end
    x, sort(idxs)
end


#Group-OMPのラッパー
function _somp!(x, y, Φ, nvar, idxs; kwargs...)
    γ = get(kwargs, :γ, 1.0) # しきい値の倍率
    max_iter = get(kwargs, :max_iter, size(y,1)) # 繰り返し回数

    if size(y,2) == size(y,3) == 1
        x = vec(x); y = vec(y)
    else
        x = reshape(x, size(x,1), :)
        y = reshape(y, size(y,1), :)
    end

    n, l = size(y,1), size(y,2)
    if length(nvar)==1
        ϵ = γ * sqrt(nvar) * sqrt(n*l + 2*sqrt(n*log(n)))
    else
        ϵ = γ * sqrt(mean(nvar)) * sqrt(n*l + 2*sqrt(n*log(n)))
        # ϵ = γ * sqrt(mean(nvar)) * sqrt(l*(n + 2*sqrt(n*log(n))))
    end
    x, idxs = Greedy.somp!(x, y, Φ, idxs; ϵ=ϵ, max_iter=max_iter)
    x, idxs
end

# SBLのラッパー
function _sbl!(x, y, Φ, nvar, idxs; args...)
    if length(nvar)==1
        x = reshape(x, size(x,1), :)
        y = reshape(y, size(y,1), :)
    else
        args = merge(args, (β = 1 ./ nvar,))
        y = reshape(y, size(y,1), :)
        x = reshape(x, size(x,1), :)
    end
    ϵ = 0.1 # 停止条件
    max_iter = get(args, :max_iter, size(y,1))
    β = 1/nvar
    model = SBL.Params(y, Φ, β=β)
    model = SBL.sbl(model, ϵ=ϵ, max_iter=max_iter)
    x = model.μ
    return x
end

# BMPのラッパー
function _bmp!(x, y, Φ, nvar, idxs; args...)
    if length(nvar)==1
        x = reshape(x, size(x,1), :)
        y = reshape(y, size(y,1), :)
    else
        y = reshape(y, size(y,1), :)
        x = reshape(x, size(x,1), :)
    end
    idxs = isnothing(idxs) ? Int[] : idxs
    ϵ = 0.0
    θ = get(args, :θ, 0.5)
    max_iter = get(args, :max_iter, size(y,1))
    β = iszero(nvar) ? 1e10 : 1/mean(nvar)
    model = BMP.Model(y, Φ, θ=θ, β=β, max_iter=max_iter)
    !isempty(idxs) && push!(model.gs, idxs)
    x = BMP.bmp!(x; model=model, ϵ=ϵ)
    idxs = union!(idxs, sort!(model.gs[1]))
    # @show model.gs
    x, idxs
end


function _interp_ls!(x, y, nvar; args...)
    tmp = @views similar(y[:,1,1]) # ベクトルコピ-
    for q = axes(y,2)
        for p = axes(y,3)
            tmp .= @views y[:,q,p]
            tmp = ifft!(tmp)
            x[:,q,p] .= tmp[1:size(x,1)]
        end
    end
    x
end

# 雑音電力密度からしきい値を決定する
function _define_epsilon(nvar::Number, y)
    iszero(nvar) && return 1e-3 # デフォルトしきい値
    N = size(y,1) # samples
    M = size(y,2)*size(y,3) # num of datasets
    return sqrt(nvar*N) * M
end
# 多次元の場合
function _define_epsilon(nvars::AbstractVector{T}, y) where T
    @assert length(nvar) == size(y,3)
    n_obs_size = size(y,1)
    n_obs_data = size(y,2)
    return sum(sqrt(nvar*n_obs_size*n_obs_data) for nvar in nvars)
end




"""
    ls(y::Vector{T}, Phi::Matrix{T}) where T

LS(最小二乗)推定
"""
function ls!(x, y, Φ)
    if size(Φ,1)==size(Φ,2)
        x[:,:] .= Φ \ y
    else
        x[:,:] .= pinv(Φ) * y
    end
    x
end
function ls!(x, y, Φ, idxs)
    if size(Φ,1) == size(Φ,2)
        x[idxs,:] .= Φ[:,idxs] \ y
    else
        x[idxs,:] .= pinv(Φ[:,idxs]) * y
    end
    x
end

"""
    lmmse(y::Vector{T}, Phi::Matrix{T}, cov::Matrix{T}) where T
    lmmse(y::Vector{T}, Phi::Matrix{T}, Ryx::Matrix{T}, Ryy::Matrix{T}, ) where T

Linear-MMSE推定
"""
function lmmse!(x, y, Φ::AbstractMatrix, nvar, idxs::Nothing, Σ::Nothing)
    if nvar isa Number
        x[:,:] .= (Φ'*Φ+ nvar*I) \ (Φ'*y)
    else
        x[:,:] .= Φ' * inv(Φ*Φ' + Diagonal(nvar)) * y
    end
    x
end
function lmmse!(x, y, Φ::AbstractMatrix, nvar, idxs::AbstractVector, Σ::Nothing)
    Φp = Φ[:,idxs]
    x[idxs,:] .= (Φp'*Φp + nvar*I) \ (Φp'*y)
    x
end
function lmmse!(x, y, Φ::AbstractMatrix, nvar, idxs::Nothing, Σ::AbstractMatrix)
    if nvar isa Number
        x[:,:] .= (Φ'*Φ+ nvar*I) \ (Φ'*y)
    else
        x[:,:] .= Φ' * inv(Φ*Φ' + Diagonal(nvar)) * y
    end
end
function lmmse!(x, y, Φ::AbstractMatrix, nvar, idxs::AbstractVector, Σ::AbstractMatrix)
    Φ_ = Φ[:,idxs]
    if nvar isa Number
        x[idxs,:] .= (Σ*Φ_'*Φ_+ nvar*I) \ (Σ*Φ'*y)
    else
        x[idxs,:] .= Φ' * inv(Φ*Σ*Φ' + Diagonal(nvar)) * y
    end
end

function lmmse_with_prior!(x, y, Rxx, Ryy)
    x .= (Rxy * inv(Ryy)) \ y
end


# ルーティン用関数
function _routine!(f, x, y, Φ, args...)
    for j in axes(y,2)
        for i in axes(x,2)
            f(x, y, Φ, args...)
        end
    end
    x
end
