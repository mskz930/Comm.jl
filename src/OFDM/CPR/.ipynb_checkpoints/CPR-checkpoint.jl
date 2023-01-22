module CPR

using FFTW, LinearAlgebra

# 矩形パルス窓関数
function u(n, g)
    if (n-g) < 0
        return 0.0
    else
        return 1.0
    end
end

# チャネル推定
function get_cfr(cir::AbstractArray{T}, nfft, cpsize) where T
    L = size(cir, 1) # チャネル次数
    cir = [cir; zeros(T, nfft-L, size(cir)[2:end]...)]
    cfr = fft!(cir, 1) # FFT
    # 補正
    if cpsize < L-1
        for q in axes(cir,2)
            for p in axes(cir,3)
                for k in 0:nfft-1
                    for l in cpsize+1:L-1
                        cfr[k+1,q,p] -= l/nfft*cir[l+1,q,p]*exp(-im*2*pi*k*l/nfft)
                    end
                end
            end
        end
    end
    return cfr
end




# CPR: 巡回修復(Cyclic Property Reconstruction)
function cpr!(r, h, X0, X1, DFTmtx, nfft, cpsize, rind, dind, D=size(r,1), domain=:freq)
    """
    arguments:　
        r, h      : 受信信号(OFDMシンボル), インパルス応答
        X1, X0    : 前判定シンボル, 現判定シンボル
        D, indices: 近似パラメータ(1~nfft), 等化するサブキャリアインデックス
    """
    if domain == :freq
        # 周波数領域で実行
        if !isnothing(X1)
            cpr_in_freq!(r, h, X1, X0, DFTmtx, nfft, cpsize, rind, dind, D)
        elseif !isnothing(X0)
            isi_remove!(r, h, X0, DFTmtx, nfft, cpsize, rind, dind, D)
        end
    elseif domain == :time
        # 時間領域で実行
        cpr_in_time!(r, h, X1, X0, F, nfft, cpsize)

    end
    return r
end

# 巡回再構築(CPR: cyclic-property-econstruction)
# 周波数領域でCPRを実行
function cpr_in_freq!(r, h, X1, X0, F, nfft, cpsize, rind, dind, D)
    maxk = maximum(dind); mink = minimum(dind);
    L = size(h,1) # チャネル次数
    c::ComplexF64 = 0
    temp::ComplexF64 = 0
    @inbounds for q in axes(r,2) # 受信アンテナ
        for k in rind # 周波数インデックス
            minr = k-D < mink ? mink : k-D
            maxr = k+D > maxk ? maxk : k+D
            for p in axes(X1,2) # 送信アンテナ
                for m in minr:maxr # 隣接キャリア
                    c=0; temp=0;
                    d = k-m # 周波数差
                    for l in cpsize+1:L-1 # 遅延軸
                        if abs2(h[l+1,p,q]) > 0
                            if k !== m
                                # ICI除去
                                for n in 0:(l-cpsize)-1 # 時間軸
                                    c += d > 0 ? F[n+1,d+1] : conj(F[n+1,-d+1])
                                    # c += exp(-im*(2*pi*d*n)/nfft)
                                end
                                temp += c * h[l+1,q,p] * F[m,l+1]
                            else
                                temp = l * h[l+1,q,p] * F[k,l+1]
                            end
                        end
                    end
                    r[k,q] += temp * X1[m,p] / nfft # 補正項を足す
                    if !isnothing(X0)
                        # ISI除去
                        r[k,q] -= temp * X0[m,p] / nfft
                    end
                end
            end
        end
    end
end

# ISI除去(周波数領域で実行)
function isi_remove!(r, h, X0, F, nfft, cpsize, rind, dind, D)
    maxk = maximum(dind); mink = minimum(dind);
    L = size(h,1) # チャネル次数
    c::ComplexF64 = 0
    temp::ComplexF64 = 0
    @inbounds for q in axes(r,2) # 受信アンテナ
        for k in rind # 周波数インデックス
            minr = k-D < mink ? mink : k-D
            maxr = k+D > maxk ? maxk : k+D
            for p in axes(X0,2) # 送信アンテナ
                for m in minr:maxr # 隣接キャリア
                    c=0; temp=0;
                    d = k-m # 周波数差
                    for l in cpsize+1:L-1 # 遅延軸
                        if abs2(h[l+1,p,q]) > 0
                            # ISI除去
                            for n in 0:(l-cpsize)-1 # 時間軸
                                c += d > 0 ? F[n+1,d+1] : conj(F[n+1,-d+1])
                                # c += exp(-im*(2*pi*d*n)/nfft)
                            end
                            temp += -c * h[l+1,q,p] * F[m,l+1]
                        end
                    end
                    r[k,q] += temp * X0[m,p] / nfft # 補正項を足す
                end
            end
        end
    end
end

# 時間領域でCPRを実行
function cpr_in_time!(r, h, X1, X0, F, nfft, cpsize)
    r = ifft!(r,1) # ifft
    L = size(h,1) # 遅延次数
    x1 = !isnothing(X1) ? ifft(X1,1) : X1
    x0 = !isnothing(X0) ? ifft(X0,1) : X0
    for q in axes(r,2)
        for p in axes(X1,2)
            for l in cpsize+1:L-1
                if abs(h[l+1]) > 0
                    if isnothing(x1)
                        r[1:l] .+= h[l+1,q,p] * view(x1,nfft-l+1:nfft, p) # ICI除去
                    end
                    if isnothing(x0)
                        r[1:l] .+= h[l+1,q,p] * view(x0,nfft-l+1:nfft, p) # ISI除去
                    end
                end
            end
        end
    end
    fft!(r,1)
end

# SINR計算
function measure_sinr(nfft, range, Pn)
    sinr = zeros(Float64,length(range)) # 配列
    Ps::ComplexF64 = 0.0
    Pisi::ComplexF64 = 0.0
    Pici::ComplexF64 = 0.0

    for k in range
        for dims in 1:ndim
            for i in 1:L
                for j in 1:L
                    # 信号電力
                    Ps = (1+(i*j)/nfft^2 - (i+j)/nfft)*h[i,dims]*conj(h[j,dims]) * exp(-2π*k(i-j)/nfft)
                    if i > G+1 && j > G+1
                        # 平均ISI電力
                        Pisi += h[i,dims]*conj(h[j,dims])* min(i,j) *  exp(-im*(2π*k*(i-j)/nfft))/nfft
                    end
                    if i > G && j > G
                        # 平均ICI電力
                        Pici += h[i,dims] * conj(h[j,dims]) * exp(-im*(2π*k*(i-j)/nfft)) * (min(i,j)-(i*j)/nfft)
                    end
                    Pn = length(nvar) > 1 ? nvar[k] : nvar
                    sinr[k] = real(Ps) / (real(Pisi) + real(Pici) + Pn)
                end
            end
        end
    end
end


# ICI行列取得
function get_ici_matrix(h::AbstractVector, nfft, cpsize)
    ici_mtx = zeros(ComplexF64, nfft, nfft)
    L = length(h) # CIR長
    1/nfft * ()
    for k in 1:nfft
        for m in in 1:nfft
            if k !== m
                for n in 1:nfft
                    c = exp(im*(2π*(m-k)*n)/nfft) # 定数
                    for l in 1:L
                        ici_mtx[k,m] += h[l] * c * exp(-im*(2π*m*l)/nfft) * u(n-l) * nfft
                    end
                end
            end
        end
    end
end

# ISI行列取得
function get_isi_matrix(h::AbstractVector, nfft, cpsize)
    isi_mtx = zeros(ComplexF64, nfft, nfft)
    L = length(h)
    for k in 1:nfft
        for m in in 1:nfft
            for n in 1:nfft
                c = exp(im*(2π*(m-k)*n)/nfft) # 定数
                for l in cpsize+2:L
                    isi_mtx[k,m] += h[l] * c * exp(-im*(2π*m*l)/nfft) * u(l-n-cpsize) * nfft
                end
            end
        end
    end
end



end # end module
