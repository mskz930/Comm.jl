# cpr_receiver.jl

# 硬判定CPR
function hard_cpr(rx_frame, tx_frame, params, indices, MAXITER; N0, cir, cfr, mod)
    nblocks = size(rx_frame, 2) # OFDMシンボル数
    temp_frame = tx_frame[:,1:size(rx_frame, 2),:,:]
    eq_frame = zeros(eltype(rx_frame), size(rx_frame))
    F = fft!(Array{ComplexF64}(I, params.nfft, params.nfft), 1) # DFT行列
    buffer = Float64[] # バッファ
    sizehint!(buffer, params.nfft*nblocks*params.Nt*qam.bps) # バッファサイズ(LLR格納)
    rind = [LinearIndices(view(indices[:frame],:,j,1))[view(indices[:frame],:,j,1) .!== 0] for j in 1:nblocks]  # 等化シンボルインデックス
    dind = [LinearIndices(view(indices[:frame],:,j,1))[view(indices[:frame],:,j,1) .< 0] for j in 1:nblocks]   # データシンボルインデックス

    iter = 0
    while iter <= MAXITER
        # チャネル推定
        cir = cir[:,:,:,:]
        cfr = CPR.correct_cfr(cfr[:,:,:,:], cir, params.nfft, params.cpsize, F) # 周波数応答取得
        # cir = ifft(cir, 1) # 平均的なインパルス応答
        # cir_est = cir_estimation() # CIRの推定
        # 雑音推定
        # nvar = CPR.measure_noisevar(cir, N0) # SINRの推定

        # CPRおよびデータ判定
        previous_symbol = nothing
        current_symbol = view(temp_frame,:,1,:)
        flag0 = false; flag1 = true
        for j in axes(rx_frame,2)
            if iter==0 && j>1
                previous_symbol = view(temp_frame,:,j-1,:)
            elseif iter>0 && j==1
                current_symbol = view(temp_frame,:,j,:)
            elseif iter>0 && j>1
                previous_symbol = view(temp_frame,:,j-1,:)
                current_symbol = view(temp_frame,:,j,:)
            end
            ici_matrix = CPR.get_ici_matrix(cir, params.nfft, params.cpsize)
            CPR.cpr!(view(eq_frame,:,j,:), view(rx_frame,:,j,:), view(cir,:,j,:,:), previous_symbol, current_symbol, F, params.nfft, params.cpsize, rind[j], dind[j], flag0, flag1) # 干渉等化(CPR)
            eq_symbol = channel_equalize(view(eq_frame,:,j,:), view(cfr,:,j,:,:), view(indices[:frame],:,j,1), 0, method=:ZF) # 信号検出
            temp_symbol = slicer(eq_symbol, :qam, M) # 暫定シンボル

            if iter <= MAXITER
                temp_symbol = qammod(mod, soft_bit)
                data_mapping!(view(temp_frame,:,j,:), temp_symbol, view(indices[:frame],j,:)) # データマッピング
            else
                append!(buffer,llrs) # バッファに溜める
            end
        end
        iter += 1
    end
    return eq_symbol
end

# 軟判定CPR
function soft_cpr(rx_frame, tx_frame, params, indices, MAXITER; N0, cir=nothing, cfr=nothing, qam, permind, decparams, decoder::Function)
    # パラメータ
    nblocks = size(rx_frame,2) # ブロック数
    eq_frame = zeros(eltype(rx_frame), size(rx_frame)) # 等化受信フレーム
    temp_frame = tx_frame[:,1:size(rx_frame,2),:] # 判定シンボル配列
    nvar = zeros(Float64,size(eq_frame)) # 雑音配列
    F = fft!(Array{ComplexF64}(I, params.nfft, params.nfft), 1) # DFT行列
    buffer = Float64[] # バッファ
    sizehint!(buffer, params.nfft*nblocks*params.Nt*qam.bps) # バッファサイズ(LLR格納)
    rind = [LinearIndices(view(indices[:frame],:,j,1))[view(indices[:frame],:,j,1) .!== 0] for j in 1:nblocks]  # 等化シンボルインデックス
    dind = [LinearIndices(view(indices[:frame],:,j,1))[view(indices[:frame],:,j,1) .< 0] for j in 1:nblocks]   # データシンボルインデックス
    zero_frame = zeros(eltype(eq_frame), size(eq_frame,1),size(eq_frame,3)) # ゼロ係数

    if !isnothing(cir)
        cfr = Ofdm.CPR.correct_cfr(cfr, cir, params.nfft, params.cpsize, F) # 周波数応答取得
        if size(cir,2)==1 && size(cir,2) !== nblocks
            cfr = repeat(cfr, 1, nblocks, 1, 1)
        end
    end
    iter = 0;
    while iter <= MAXITER
        # チャネル推定
        if isnothing(cir)
            CPR.cpr!() # パイロットシンボルに対するCPR
            #　h = cir_estimation() # CIRの推定
            # 雑音推定
            # SINRの推定
        end

        # CPRおよびデータ判定
        previous_symbol = nothing # 前シンボル初期化
        current_symbol = zero_frame # 現シンボル初期化
        current_symbol = view(temp_frame, :, 1, :)
        flag0 = false; flag1 = true # flags
        for j in axes(rx_frame,2)
            if iter==0 && j>1
                flag0 = true
                previous_symbol = view(temp_frame,:,j-1,:)
            elseif iter>0 && j==1
                flag1 = true
                current_symbol = view(temp_frame,:,j,:)
            elseif iter>0 && j>1
                flag0 = true
                previous_symbol = view(temp_frame,:,j-1,:)
                current_symbol = view(temp_frame,:,j,:)
            end
            Ofdm.CPR.cpr!(view(eq_frame,:,j,:), view(rx_frame,:,j,:), view(cir,:,j,:,:), previous_symbol, current_symbol, F, params.nfft, params.cpsize, rind[j], dind[j], flag0, flag1; N0=N0, nvar=view(nvar,:,j,:)) # 干渉等化(CPR)
            if iter==0 || iter<=MAXITER # 復号用軟値ビット生成
                eq_symbol, eq_nvar = channel_equalize(view(eq_frame,:,j,:), view(cfr,:,j,:,:), view(indices[:frame], :, j, 1), view(nvar,:,j,:), method=:MMSE, outtype=:soft) # チャネル等化
                llrs = qamdemod(qam, eq_symbol, eq_nvar) # 事前LLRs
                if iter<MAXITER
                    append!(buffer,llrs) # バッファに溜める
                end
            end
            if iter==0 && j < nblocks # ISI除去用ソフトシンボル生成
                temp_symbol = qammod(qam, llrs) # ソフトシンボル生成
                # data_remapping!(params, view(temp_frame,:,j,:),  temp_symbol, view(indices[:frame], :, j, :)) # データマッピング
            end
        end
        if iter < MAXITER
            # 誤り訂正復号
            buffer[permind] = buffer # deinterleave
            Lapp = decoder(decparams, buffer, outtype=:output) # 事後LLRs
            Lapp = view(Lapp, permind) # interleave
            temp_symbol = qammod(qam, Lapp) # ソフトシンボル生成
            # data_remapping!(params, temp_frame, temp_symbol, indices[:frame])
        end
        iter += 1
    end # end while
    return eq_frame, nvar
end
