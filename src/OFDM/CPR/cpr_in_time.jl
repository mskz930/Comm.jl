# cpr.jl: 時間領域CPR等化


# CPRのラッパー関数
function cpr(ofdm::Ofdm, r::AbstractArray{T}, h, inds, slice; mehtod="cpr", n_iter, flag=false) where T
    r = reshape(r, ofdm.n_Ts, :, ofdm.n_dims[2])
    N, G, K = ofdm.n_fft, ofdm.n_gi, size(r,2)
    y = similar(r)                               # 等化受信信号(時間領域)
    Y = Array{T}(undef, N, K)                    # 等化受信信号(周波数領域)
    X = zeros(T, N, K)                           # 判定シンボルの格納配列
    x = Array{T}(undef, N+G, K)
    FFT  = @views plan_fft(x[G+1:end,1], 1)
    IFFT = @views plan_ifft(x[G+1:end,1], 1)

    # チャネル推定
    H  =  OFDM.ChEst.convert_to_cfr(h, N, G)

    # 分岐
    if method == "CPR"
        # channel_compensation!(H, h, N, G)            # チャネル係数補正
        _cpr!(x, X, y, Y, h, H, r, FFT, IFFT, inds, slice, N, G, K, n_iter)
    elseif method == "W-CPR"
        w = L-G>0 ? zeros(L-G) : zeros(0)
        _weight_calc!(w, h, G)
        _weighted_cpr!(x, X, y, Y, h, H, r, FFT, IFFT, inds, slice, N, G, K, n_iter, flag)
    elseif method == "SC-CPR"

    end
    X
end

# 時間領域CPRの実行関数
function _cpr!(x, X, y, Y, h, H, r, FFT, IFFT, inds, slice, N, G, K, n_iter, flag)
    L = size(h,1)-1
    iter = 0
    while iter < n_iter
        y .= r
        for k in 1:K
            # ISI subtraction
            if k > 1
                x[(G+1):end,k-1] = @views (IFFT * X[:,k-1]) .* sqrt(N)
                @views _isi_remove!(y[:,k], x[:,k-1], h, G, L)
            end

            # Cyclic Reconstruction
            if (iter > 0) && k < size(r,2)
                x[G+1:end,k] .= (IFFT * @view X[:,k]) * sqrt(N)  # 推定OFDM時間シンボル
                @views _ici_remove!(y[:,k], x[:,k], h, G, L)
            end

            # Data detection
            Y[:,k] = (FFT * @view y[G+1:end,k]) / sqrt(N)
            @views _data_detection!(X[:,k], Y[:,k], H, inds[:,k], slice)
        end
        iter += 1
    end
end


# weighted_cprの実行部分
function _weighted_cpr!(x, X, y, Y, h, H, r, w, FFT, IFFT, inds, slice, N, G, K, n_iter, flag)
    L = size(h,1)-1
    iter = 0
    while iter < n_iter
        y .= r
        for k in 1:K
            # isi subtraction
            if k > 1
                x[:,k-1] = @views (IFFT * X[:,k-1]) .* sqrt(N)
                @views _isi_remove!(y[:,k], x[:,k-1], h, G, L)
            end

            # Weighted sum / Cyclic Reconstruction
            if (iter == 0) && k < size(r,2) && flag
                @views _weighted_sum!(y[:,k], r[:,k+1], w, G, L)

            elseif (iter > 0) && k < size(r,2)
                x[:,k] .= (IFFT * @view X[:,k]) .* sqrt(N)  # 推定OFDM時間シンボル
                @views _ici_remove!(y[:,k], x[:,k], h, G, L)
            end

            # data_detection
            Y[:,k] = (FFT * @view y[G+1:end,k]) / sqrt(N)
            @views _data_detection!(X[:,k], Y[:,k], H, inds[:,k], slice)
        end
        iter += 1
    end
end

# 重みベクトルの計算
function _weight_calc!(w, h, G)
    L = size(h,1)-1
    for i in eachindex(w)
        for l in G+i:L
            w[i] += abs2(h[l+1])
        end
    end
    w ./= sum(abs2(h[i]) for i in eachindex(h))
end

# 重み付きICI補償
function _weighted_sum!(y, r, w, G, L)
    for i in 1:L-G
        y[G+i] += w[i] * r[i]
    end
end

# 巡回シフト
function _sccpr!()
    L = size(h,1)-1
    iter = 0
    while iter < n_iter
        y .= r
        for k in 1:K
            # ISI subtraction
            if k > 1
                x[(G+1):end,k-1] = @views (IFFT * X[:,k-1]) .* sqrt(N)
                @views _isi_remove!(y[:,k], x[:,k-1], h, G, L)
            end

            # Cyclic Reconstruction
            if flag && iter==0
                @views _sample_estimate!(x[:,k], y[:,k], h, N, G, L)
                @views _ici_remove!(y[:,k], x[:,k], h, G, L)

            elseif iter > 0 && k < size(r,2)
                x[G+1:end,k] .= (IFFT * @view X[:,k]) * sqrt(N)  # 推定OFDM時間シンボル
                @views _ici_remove!(y[:,k], x[:,k], h, G, L)
            end

            # Data detection
            Y[:,k] = (FFT * @view y[G+1:end,k]) / sqrt(N)
            @views _data_detection!(X[:,k], Y[:,k], H, inds[:,k], slice)
        end
        iter += 1
    end
end


# 逐次サンプル推定
function _sample_estimate!(x, r, h, N, G, L)
    nzinds = findall(!iszero, vec(h))
    fst = findfirst(!iszero, h)
    drops = 1

    for i in 1:N+G
        x[i] = r[i]
        if i < L+1
            for j in Iterators.drop(nzinds,drops)
                j>i  && break
                x[i] -= x[i-(j-1)] * h[j]
            end
        else
            for j in Iterators.drop(nzinds,drops)
                x[i] -= x[i-(j-1)] * h[j]
            end
        end
        x[i] /= h[fst]
        # @show i, x[i]
    end
end

function _time_domain_filter_cpr!()
    iter = 0
    while iter < n_iter
        y .= r
        for k in 1:K
            # ISI subtraction
            if k > 1
                x[(G+1):end,k-1] = @views (IFFT * X[:,k-1]) .* sqrt(N)
                @views _isi_remove!(y[:,k], x[:,k-1], h, G, L)
            end

            # Cyclic Reconstruction
            if flag && iter==0
                
                @views _ici_remove!(y[:,k], x[:,k], h, G, L)

            elseif iter > 0 && k < size(r,2)
                x[G+1:end,k] .= (IFFT * @view X[:,k]) * sqrt(N)  # 推定OFDM時間シンボル
                @views _ici_remove!(y[:,k], x[:,k], h, G, L)
            end

            # Data detection
            Y[:,k] = (FFT * @view y[G+1:end,k]) / sqrt(N)
            @views _data_detection!(X[:,k], Y[:,k], H, inds[:,k], slice)
        end
        iter += 1
    end
end



# ISI除去
function _isi_remove!(y, x̂, ĥ, G, L)
    for l in 1:L
        if !iszero(ĥ[l+1])
            for i in 1:l
                y[i] -= ĥ[l+1] * x̂[end-l+i]
            end
        end
    end
end

# ICI除去
function _ici_remove!(y, x̂, ĥ, G, L)
    for l in G+1:L
        if !iszero(ĥ[l+1])
            for i in 1:l-G
                y[i+G] += ĥ[l+1] * x̂[end-l+i]
            end
        end
    end
end

# データ検出
function _data_detection!(X, Y, H, inds, slice)
    for i in eachindex(inds)
        if inds[i]<0
            X[i] = slice(Y[i] / H[i])
        end
    end
end
