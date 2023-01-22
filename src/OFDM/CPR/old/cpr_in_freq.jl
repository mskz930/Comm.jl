# cpr_in_freq.jl

function cpr_in_freq()
    
end

"""
ハードシンボルによるCPRの実行
"""
function cpr_in_freq!(Y, Hi, X0)
    # ISI除去
    flag0 && _isi_remove_in_freq!(Y, Hi, X0, inds0, D)

    # ICI除去
    flag1 && _ici_remove_in_freq!(Y, Hi, X1, inds1, D)
end

# ISI除去
function _isi_remove_in_freq!(Y, Hi, X0, inds, D)
    N = size(Y,1)
    st = en = 0
    for q in axes(Y,2) # Nr
        for p in axes(X0,2) # Nt
            for k in axes(Y,1) # freq
                st = ifelse(k-D<0, 1, k-D) # k-D (if k-D>1)
                en = ifelse(k-D>N, N, k+D) # k-D (if k+D<N)
                for m in st:en # freq
                    if inds[k,p] != 0
                        Y[k,q] -= Hi[k,m,q,p] * X0[m,p]
                    end
                end
            end
        end
    end
end

# ICI除去(補正)
function _ici_remove_in_freq!(Y, Hi, X1, inds, D)
    N = size(Y,1)
    st = en = 0
    for q in axes(Y,2) # Nr
        for p in axes(X0,2) # Nt
            for k in axes(Y,1) # freq
                for m in st:en # freq
                    if k != m && inds[k,p] != 0
                        Y[k,q] += Hi[k,m,q,p] * X1[m,p]
                    end
                end
            end
        end
    end
end




"""
ソフトシンボルによるCPRの実行
"""
function soft_cpr_in_freq!(Y, nvar, Hi, X0, X1, inds0, inds1, D=size(Y,1), flag0, flag1)
    # ISI除去
    flag0 && _soft_isi_remove!(Y, nvar, Hi, X0, inds0, D)

    # ICI除去
    flag1 && _soft_ici_remove!(Y, nvar, Hi, X1, inds1, D)

    # 等価雑音電力のみ計算
    (!flag0&&!flag2) && _nvar_calc(Y, nvar, Hi, inds1)
end

function _soft_isi_remove_in_freq!(Y, nvar, Hi, X0, inds, D)
    N = size(Y,1)
    st = en = 0
    for q in axes(Y,2) # Nr
        for p in axes(X0,2) # Nt
            for k in axes(Y,1) # freq
                st = ifelse(k-D<0, 1, k-D) # k-D (if k-D>1)
                en = ifelse(k-D>N, N, k+D) # k-D (if k+D<N)
                for m in st:en # freq
                    if inds[k,p] != 0
                        Y[k,q] -= Hi[k,m,q,p] * X0[m,p].mean
                        nvar[k,q] -= abs2(Hisi[k,m,q,p]) * X0[m,p].var
                    end
                end
            end
        end
    end
end

function _soft_ici_remove_in_freq!(Y, var, Hi, X1, D, inds)
    N = size(Y,1)
    st = en = 0
    for q in axes(Y,2) # Nr
        for p in axes(X0,2) # Nt
            for k in axes(Y,1) # freq
                for m in st:en # freq
                    if k != m && inds[k,p] != 0
                        Y[k,q] += Hi[k,m,q,p] * X1[m,p].mean
                        nvar[k,q] += abs2(Hisi[k,m,q,p]) * X1[m,p].var
                    end
                end
            end
        end
    end
end

# 干渉＋雑音電力のみ算出
function _nvar_calc_in_freq!(Y, nvar, Hi, inds)
    for q in axes(Y,2) # Nr
        for p in axes(X,2) # Nt
            for k in axes(Y,1) # freq
                for m in axes(Y,1) # freq
                    if inds[k,p] != 0 && k != m
                        nvar[k,q] += abs2(Hisi[k,m,q,p])
                    end
                end
            end
        end
    end
end




"""
等価雑音分散の計算
"""
function equiv_noisevar_calc!(nvar, Hici=nothing, Hisi=nothing, sinds=nothing, inds1=nothing, inds0=nothing, D=size(nvar,1))
    fft_size, Nt, Nr = size(Hici,1), size(Hici,3), size(Hici,4)
    for q in 1:Nr
        for k in sinds # frequency index in current symbol frame
            ivar = 0.0 # 干渉電力
            for p in 1:Nt
                #  isi power
                if !isnothing(Hisi) && !isnothing(inds0)
                    for m in inds0
                        if inds0[m,p] !== 0
                            ivar += abs2(Hisi[k,m,q,p]) # ICI残差電力
                        end
                    end
                end
                # ici power
                if !isnothing(Hici) && !isnothing(inds1)
                    for m in 1:fft_size
                        if (m!==k) && inds1[m,p] !== 0
                            (ivar += abs2(Hici[k,m,q,p])) # ISI残差電力
                        end
                    end
                end
            end
            nvar[k,q] = ivar
        end
    end
end



"""

巡回修復(Cyclic Property Reconstruction)
"""
function cpr!(eq, r, h, Xp, Xc, DFTmtx, nfft, cpsize, rind, dind, flag0, flag1; D=size(r,1), nvar=nothing, N0=nothing, domain=:freq)
    """
    arguments:　
        r, h      : 受信信号(OFDMシンボル), インパルス応答
        Xc, Xp    : 前判定シンボル, 現判定シンボル
        D, indices: 近似パラメータ(1~nfft), 等化するサブキャリアインデックス
    """
    if domain == :freq # 周波数領域で実行
        cpr_in_freq!(eq, r, h, Xc, Xp, N0, nvar, DFTmtx, nfft, cpsize, rind, dind, flag0, flag1, D)
    elseif domain == :time # 時間領域で実行
        cpr_in_time!(eq, r, h, Xc, Xp, nfft, cpsize)
    end
end

# 巡回再構築(CPR: cyclic-property-econstruction)
# cpr_in_freq!: 周波数領域でICI除去・巡回再構築を実行
function cpr_in_freq!(eq, r, h, Xc, Xp, N0, nvar, F, nfft, cpsize, rind, dind, flag0, flag1, D)
    maxk = maximum(dind) # データインデックスの最大値
    mink = minimum(dind) # データインデックスの最小値
    L = size(h,1) # チャネル次数
    pind = findall(x -> x > 0, abs.(view(h,:,1,1,1)))
    pind = filter(x -> x > cpsize+1, pind) .- 1 # 非ゼロなパスインデックス
    @inbounds for q in axes(r,2) # 受信アンテナ
        for k in rind # 周波数インデックス
            minr = k-D < mink ? mink : k-D
            maxr = k+D > maxk ? maxk : k+D
            val = zero(ComplexF64)
            ipow = zero(Float64)
            for p in axes(h,3)
                for m in minr:maxr # 隣接キャリア
                    coef = zero(ComplexF64) # チャネル係数
                    if k !== m
                        d = k-m # 周波数差
                        for l in pind # 遅延軸
                            b = zero(ComplexF64) # 位相
                            if d >= 0
                                for n in 0:(l-cpsize)-1 # 時間軸
                                    b += F[n+1,d+1]
                                end
                            else
                                for n in 0:(l-cpsize)-1 # 時間軸
                                    b += conj(F[n+1,-d+1])
                                end
                                # b += exp(-im*2*pi*d*n/nfft)
                            end
                            coef += b * h[l+1,q,p] * F[m,l+1] # ICI係数
                            # coef += b * h[l+1,q,p] * exp(-im*2*pi*(m-1)*l/nfft)
                        end
                    else
                        for l in pind # 遅延軸
                            coef += l * h[l+1,q,p] * F[k,l+1]# 巡回損失補正係数
                        end
                    end
                    if flag1 && m!==k
                        val -= coef * Xc[m,p] / nfft # ICI成分
                    end
                    if flag0
                        val += coef * Xp[m,p] / nfft # ISI成分
                    end
                    if !isnothing(nvar) && !isnothing(Xc) && m!==k
                        ipow += real(abs2(coef) * (1 - abs2(Xc[m,p])) / nfft^2) # ICI残差電力
                    end
                    if !isnothing(nvar) && !isnothing(Xp)
                        ipow += real(abs2(coef) * (1 - abs2(Xp[m,p])) / nfft^2) # ISI残差電力
                    end
                end
                eq[k,q] = r[k,q] - val
                if !isnothing(nvar)
                    nvar[k,q] = N0 + ipow
                end
            end
        end
    end
end

function soft_cpr!(Y, X1, X0, Hici, Hisi, N=nothing; yinds, x0inds, x1inds, D=size(Y,1))
    if !isnothing(X1) # ICI除去
        for q in axes(Y,2)
            for p in axes(X1,2)
                @views _remove!(Y, X1, Hici[:,:,q,p], yinds, x1inds, D)
            end
        end
    end
    if !isnothing(X1) # ISI除去
        for q in axes(Y,2)
            for p in axes(X1,2)
                @views _remove!(Y, X0, Hisi[:,:,q,p], yinds, x0inds, D)
            end
        end
    end
end

# 干渉除去
function _remove!(Y, X, Hi, W, D=size(y,1); yinds, xinds)
    #  Y: 受信ベクトル
    #  X: 送信ベクトル(推定)
    # Hi: 干渉係数行列
    #  W: 雑音電力ベクトル
    #  D: 探索範囲制限定数
    for m in yinds
        val  =  complex(0); pow = float(0)
        for k in xinds
            if m-D > k || k > m+D
                continue
            else
                val += Hi[m,k] * X[k].mean
                pow  += X[m].var * abs2(Hi[k,m,q,p])
            end
        end
        Y[m] -= val;
        W[m] += pow;
    end
end
