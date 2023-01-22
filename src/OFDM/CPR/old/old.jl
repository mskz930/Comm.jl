
"""
干渉電力計算
"""
function ipow_calc!(Xp, h, nvar, F, nfft, cpsize, rind, dind, D)
    maxk = maximum(dind); mink = minimum(dind);
    L = size(h,1) # チャネル次数
    c::ComplexF64 = 0
    coef::ComplexF64 = 0
    @inbounds for q in axes(r,2) # 受信アンテナ
        for k in rind # 周波数インデックス
            minr = k-D < mink ? mink : k-D
            maxr = k+D > maxk ? maxk : k+D
            val = zero(eltype(eq))
            for p in axes(Xp,2) # 送信アンテナ
                for m in minr:maxr # 隣接キャリア
                    c=0; coef=0;
                    d = k-m # 周波数差
                    for l in cpsize+1:L-1 # 遅延軸
                        if abs2(h[l+1,p,q]) > 0
                            # ISI除去
                            for n in 0:(l-cpsize)-1 # 時間軸
                                c += d > 0 ? F[n+1,d+1] : conj(F[n+1,-d+1])
                                # c += exp(-im*(2*pi*d*n)/nfft)
                            end
                            coef += -c * h[l+1,q,p] * F[m,l+1]
                        end
                    end
                    val -= coef * Xp[m,p] / nfft # ISI除去
                    if !isnothing(nvar)
                        nvar[k,q] += abs2(coef) * (1 - abs2(Xp[m,p])) / nfft^2 # ISI項
                    end
                end
            end
            eq[k,q] = r[k,q] - val
        end
    end
end



# measure_noise: 等化雑音推定
function measure_noise(h, nfft, cpsize, Nr, N0)
    L = size(h, 1) # チャネル次数
    N = nfft
    G = cpsize
    Pisi = zeros(nfft, size(h,2), Nr)
    Pici = zeros(nfft, size(h,2), Nr)
    for j in 1:size(h,2)
        for k in 1:Nr
            if j==1
                Pisi[:,j,k], Pici[:,j,k] = ipow_calc(h, nfft, cpsize, false)
            else
                Pisi[:,j,k], Pici[:,j,k] = ipow_calc(h, nfft, cpsize, true)
            end
        end
    end
    Pisi, Pici
end

# get_ipow: 干渉電力計算
function get_ipow(h::AbstractVector{T}, N, G, flag) where T
    L = length(h)
    Pici = zeros(Float64,N)
    Pisi = zeros(Float64,N)
    for k in 1:N
        for p in 2:L
            for q in 2:L
                if h[p]!==zero(T) && h[q]!==zero(T)
                    c = exp(-im*2*pi*(p-q)*k/N)
                    Pici[k] += real(h[p]*conj(h[q])*c*(min(p,q)/N-(p*q/N^2))) # ICI
                    if flag && (p>G+1) && (q>G+1)
                        Pisi[k] += real(h[p]*conj(h[q])*c*min(p,q)/N) # ISI
                    end
                end
            end
        end
    end
    return Pisi, Pici
end

"""
    measure_sinr()
SINR計算
"""
function measure_sinr(nfft, range, N0)
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
