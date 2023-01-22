# cpr_receiver.jl

# 硬判定CPR
function hard_cpr!(rx_frame, tx_frame, cir, N0, MAXITER, ofdmparams, indices, mod)
    # チャネル推定
    cfr = Ofdm.CPR.get_cfr(cir, params.nfft, params.cpsize) # 周波数応答取得
    cir = ifft(cir, 1) # 平均的なインパルス応答
    # cir_est = cir_estimation() # CIRの推定
    # 雑音推定
    # nvar = CPR.measure_noisevar(cir, N0) # SINRの推定

    # CPRおよびデータ判定
    previous_symbol = nothing; current_symbol = nothing
    for j in axes(rx_frame,2)
        if iter==0 && j>0
            previous_symbol = view(temp_frame,:,j-1,:)
        elseif iter>0 && j==0
            current_symbol = view(temp_frame,:,j,:)
        elseif iter>0 && j>0
            previous_symbol = view(temp_frame,:,j-1,:)
            current_symbol = view(temp_frame,:,j,:)
        end
        ici_matrix = CPR.get_ici_matrix(cir, params.nfft, params.cpsize)
        cpr!(view(eq_frame,:,j,:), cir, previous_symbol, current_symbol, rind, dind) # 干渉等化(CPR)
        eq_symbol = channel_equalize(view(eq_frame,:,j,:), view(cfr,:,j,:,:), frameinfo, view(nvar,:,j,:), method=moethod) # 信号検出
        if iter < MAXITER
            temp_symbol = slicer(eq_symbol, :qam, M) # 暫定シンボル
            temp_symbol = qammod(mod, soft_bit)
            data_mapping!(view(temp_frame,:,j,:), temp_symbol, view(frameinfo,j,:)) # データマッピング
        else
            append!(buffer, eq_symbol) # バッファに溜める
        end
    end
    return bufferd
end

# 軟判定CPR
function soft_cpr(rx_frame, tx_frame, params, frameinfo, MAXITER; N0, w_size=size(rx_frame,2), cir=nothing, cfr=nothing, qam, permind, decparams, decoder)

    # 初期パラメータ
    nblocks = size(rx_frame,2) # ブロック数
    eq_frame = zeros(eltype(rx_frame), size(rx_frame)) # 等化受信フレーム
    temp_frame = tx_frame[:,1:size(rx_frame,2),:] # 判定シンボル配列
    zero_frame = zeros(eltype(eq_frame), size(eq_frame,1),size(eq_frame,3)) # ゼロ係数
    nvar = zeros(Float64,size(eq_frame)) # 雑音配列
    F = fft!(Array{ComplexF64}(I, params.nfft, params.nfft), 1) # DFT行列
    rind = LinearIndices(view(frameinfo,:,1,1))[view(frameinfo,:,1,1) .!== 0] # 等化シンボルインデックス
    dind = [LinearIndices(view(frameinfo,:,j,1))[view(frameinfo,:,j,1) .< 0] for j in 1:nblocks] # データシンボルインデックス
    pind = [LinearIndices(view(frameinfo,:,j,1))[view(frameinfo,:,j,1) .> 0] for j in 1:nblocks] # パイロットインデックス

    # 周波数応答計算
    if !isnothing(cir)
        cfrc = correct_cfr(copy(cfr), cir, params.nfft, params.cpsize, F) # 周波数応答取得
        if size(cir,2)==1 && size(cir,2) !== nblocks
            cfrc = repeat(cfr_, 1, nblocks, 1, 1)
        end
        cfrc_ = zeros(eltype(cfrc),size(cfrc))
    else
        # channel_estimate() # 初期チャネル推定
        # N0 # 雑音推定
    end

    # 繰り返し計算
    iter = 0;
    while iter <= MAXITER
        # 出力バッファ
        buffer = Float64[]
        sizehint!(buffer, params.nfft*nblocks*params.ndim[1]*qam.bps) # バッファサイズ(LLR格納)

        # チャネル推定
        if iter>0
        end

        # CPR&データ検出
        previous_symbol = nothing # 前シンボル列
        current_symbol = zero_frame # 現シンボル列
        flag0 = false; flag1 = true # flags
        div = Int(size(rx_frame,2)/w_size) # フレーム分割数
        for p in 1:div
            ici_mtx, isi_mtx = CPR.get_coef_matrix(view(cir,:,1,:,:), params.nfft, params.cpsize, F) # 干渉行列生成
            for j in Int(1+(p-1)*w_size):Int(p*w_size)
                current_symbol = view(temp_frame,:,j,:)
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
                # cpr!(view(eq_frame,:,j,:), view(rx_frame,:,j,:), view(cir,:,j,:,:), previous_symbol, current_symbol, F, params.nfft, params.cpsize, rind[j], dind[j], flag0, flag1; D=16, N0=N0, nvar=view(nvar,:,j,:), domain=:freq) # 干渉等化(CPR)
                cpr!(view(eq_frame,:,j,:), view(rx_frame,:,j,:), ici_mtx, isi_mtx, previous_symbol, current_symbol, flag0, flag1, rind, rind, nvar=view(nvar,:,j,:), N0=N0)

                if iter==0 || iter<=MAXITER # 復号用軟値ビット生成
                    if iter==0
                        eq_symbol, eq_nvar = channel_equalize(view(eq_frame,:,j,:), view(cfrc,:,j,:,:), view(frameinfo, :, j, 1), view(nvar,:,j,:), method=:ZF, outtype=:soft) # チャネル等化
                    else
                        eq_symbol, eq_nvar = channel_equalize(view(eq_frame,:,j,:), view(cfr,:,j,:,:), view(frameinfo, :, j, 1), view(nvar,:,j,:), method=:ZF, outtype=:soft) # チャネル等化
                    end
                    llrs = qamdemod(qam, eq_symbol, eq_nvar) # 事前LLRs
                    if iter<=MAXITER
                        append!(buffer,llrs) # バッファに溜める
                    end
                end
                if iter==0 && j<=nblocks # ソフトシンボル生成
                    # temp_symbol = qammod(qam, llrs)
                    # data_remapping!(params, view(temp_frame,:,j,:),  temp_symbol, view(frameinfo, :, j, :)) # データマッピング
                end
            end
        end
        if iter < MAXITER
            if isnothing(decoder)
                temp_symbol = qammod(qam, buffer) # ソフトシンボル生成
                data_remapping!(params, temp_frame, temp_symbol, frameinfo) # シンボルマッピング
            else
                buffer[permind] = buffer # deinterleave
                Lapp = decoder(decparams, buffer, outtype=:output) # 事後LLRs
                Lapp = Lapp[permind] # interleave
                temp_symbol = qammod(qam, Lapp) # ソフトシンボル生成
                data_remapping!(params, temp_frame, temp_symbol, frameinfo) # シンボルマッピング
            end
        end
        iter += 1
    end # end while
    return eq_frame, nvar
end

# チャネル推定
#=
function channel_estimation()
    ave_cfr = get_ave_cfr(cfr) # チャネル応答の平均化
    N0 = noise_est() # ノイズ推定
    tap_indices = nonzero_tap_detection() # 非ゼロタップ検出
    cfr_est =  # LMMSE推定
end

# チャネル再推定
function channel_re_estimation()
    cpr
end
=#


function hard_cpr(params, rx_frame, MAXITER)
    eq_frame = zeros(eltype(rx_frame), size(rx_frame)) # 等化受信フレーム
    temp_frame = tx_frame[:,1:size(rx_frame,2),:] # 判定シンボル配列
    nvar = zeros(Float64,size(eq_frame)) # 雑音配列
    F = fft!(Array{ComplexF64}(I, params.nfft, params.nfft), 1) # DFT行列
    rind = LinearIndices(view(frameinfo,:,j,1))[view(frameinfo,:,1,1) .!== 0]  # シンボルインデックス
    dind = [LinearIndices(view(frameinfo,:,j,1))[view(frameinfo,:,j,1) .< 0] for j in 1:nblocks]  # データシンボルインデックス
    pind = [LinearIndices(view(frameinfo,:,j,1))[view(frameinfo,:,j,1) .> 0] for j in 1:nblocks] # パイロットインデックス
    iter=0
    while iter < MAXITER
        hard_cpr!(params, eq_frame, rx_frame, temp_frame, nvar, N0, mod, iter, indices, F, rind, dind, pind) # softCPR

        iter += 1
    end
    buffer
end
