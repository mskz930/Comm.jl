
# soft_cpr関数の実行部分
function soft_cpr!(Y, X, nvar, CFR, Hisi, Hici, ofdm, modtype, dinds, pinds, iter; D=ofdm.nfft, N0, CIR=nothing, chest=(flag=false, ), detection=:MMSE) # , flag_chest, chest, detection, iter)
    buffer = Float64[]
    # (iter==0 && chest[:flag])

    if chest[:flag] && iter > 0 # パイロットシンボルに対するCPR
        pilot_time_inds = gen_time_indices(ofdm.pilot, size(Y,2)) # パイロットシンボルの送信タイムテーブル
        for n in axes(Y,2)
            pilot_time_inds[n]==0 && continue
            if n > 1
                @views _soft_isi_remove!(Y[:,n,:], nvar[:,n,:], Hisi, X[:,n-1,:], pinds[n], D)
            end
            @views _soft_ici_comp!(Y[:,n,:], nvar[:,n,:], Hici, X[:,n,:], pinds[n], D)
        end
        # チャネル再推定
        ave_pilots = Chest.pilot_symbol_average(ofdm, Y)
        CIR = Chest.cir_estimate(ofdm, ave_pilots, N0=N0, method=chest[:method], L=chest[:L], args=chest[:args])
        Hisi[:] .= zero(eltype(Hisi)); Hici[:] .= zero(eltype(Hici))
        gen_coef_matrix!(Hisi, Hici, CIR, ofdm.nfft, ofdm.ngi, ofdm.ndims, domain=:freq)
    end

    # データシンボルに対するCPR
    for n in 1:size(Y,2)
        # ICI, ISI除去
        n > 1 && @views _soft_isi_remove!(Y[:,n,:], nvar[:,n,:], Hisi, X[:,n-1,:], dinds[n], D)
        @views _soft_ici_comp!(Y[:,n,:], nvar[:,n,:], Hici, X[:,n,:], dinds[n], D)

        # チャネル等化
        eq_data, eq_nvar = @views channel_equalize(Y[dinds[n],n,:],
                                                   CFR[dinds[n],:,:],
                                                   n_dims=ofdm.ndims,
                                                   N0=nvar[dinds[n],n,:],
                                                   method=:MMSE,
                                                   output=:soft)

        # 軟判定
        demod_LLRs = demod(modtype, eq_data, eq_nvar)
        append!(buffer, vec(demod_LLRs)) # buffer
        if n < size(Y,2) && iter==0
            soft_data = mod(modtype, demod_LLRs)
            @views _remap!(X[:,n,:], soft_data, dinds[n])
        end
    end
    buffer
end


# ISI除去
function _soft_isi_remove!(Y, nvar, Hisi, X, yinds, D)
    N = size(Y,1)
    fst = lst = 0
    for k in yinds # freq
        fst = ifelse(k-D<=0, 1, k-D) # k-D (if k-D>1)
        lst = ifelse(k+D>N, N, k+D) # k-D (if k+D<N)
        for m in fst:lst
            for q in axes(Y,2) # Nrx
                for p in axes(X,2) # Ntx
                    Y[k,q]    -= Hisi[k,m,q,p] * X[m,p].mean
                    nvar[k,q] += abs2(Hisi[k,m,q,p]) * X[m,p].var
                end
            end
        end
    end
end

# ICI補償
function _soft_ici_comp!(Y, nvar, Hici, Xc, yinds, D)
    N = size(Y,1)
    for k in yinds # freq
        fst = ifelse(k-D<=0, 1, k-D) # k-D (if k-D>1)
        lst = ifelse(k+D>N, N, k+D)  # k-D (if k+D<N)
        for m in fst:lst
            k == m && continue
            for q in axes(Y,2) # Nrx
                for p in axes(Xc,2) # Ntx
                    Y[k,q] += Hici[k,m,q,p] * Xc[m,p].mean
                    nvar[k,q] += abs2(Hici[k,m,q,p]) * Xc[m,p].var
                end
            end
        end
    end
end

_prod(c, X) = c * X


# シンボルマップ
function _remap!(X, data_seq, inds)
    Ntx = size(X,2)
    data_seq = reshape(data_seq, Ntx, :)

    for (idx,i) in enumerate(inds)
        for j in 1:Ntx
            X[i,j] = data_seq[j,idx]
        end
    end
end
