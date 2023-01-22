

function soft_cpr_equalizer(ofdm::Ofdm, rx_frame::OfdmFrame, )
    # 配列確保
    eq_data = similar(rx_data)
    temp_data = Array{SoftSymbol}(undef, size(tx_data))
    buffer = Float64[]; sizehit!(max_buf_size)

    # 事前計算
    chest      = channel_estimate(rx_data, tx_data, findex) # チャネル推定
    Hici, Hisi = channel_coef_mtx(chest)                    # 干渉係数行列
    data_inds  = get_data_indices(rx_frame)                 # 
    pilot_inds = get_pilot_indices(rx_frame)                #

end

# CPRの内容
function _soft_cpr_equalizer!(rx_sig, CFR, )
    iter = 0; flag=true
    while iter < n_iters && flag
        # copy
        eq_data .= r

        # Soft-Feedback
        if iter > 0
            interleave!(decode_llrs, permind)
            mod(modparms, decoded_llrs)
        end

        # パイロットに対するCPRの実行
        iter > 0 && _pilot_part()

        for n in 1:framelen
            Y  = @view rx_data[:,n,:]
            H  = @view CFRs[:,n,:,:]
            X0 = n == 0 ? nothing : @view X0[:,n,:]
            X1 = @views txdata[:,n-1,:], txdata[:,n,:]

            # 干渉除去
            soft_cpr!(Y, X1, X0, Hici, Hisi)

            # 信号検出
            demod_sig  = signal_detection(Y, H, method=method)

            # LLR判定
            demod_LLRs = demod(modtype, demod_sig)

            # BufferにLLRを保存
            append!(buffer, demod_LLRs)

            # ソフトシンボル生成->マッピング
            if n < framelen
                soft_symbols = mod(modtype, decoded_LLRs)
                # data_mapping!(X0, soft_symbols, view(findex,:,j))
            end
        end

        # SISO復号器
        if !isnothing(siso_decoder)
            deinterleave!(buffer, permind)
            decoded_llrs, flag = _siso_decoder(buffer, decoder)
        end
        iter += 1
    end
    decoded_LLRs
end

#
function _data_part()

end

#
function _pilot_part()
    for n in axes(rx_sig, 2)
        if pilot_inds["time"][n] > 0
            dinds = data_inds[n]
            pinds = pilot_inds["freq"][n]
            soft_cpr!(Y, X0, X1, yinds=pinds, x0inds=dinds[n-1], x1inds=dinds[n])
        end
    end
    chest = channel_restimate()
    Hici, Hisi = coef_mtx_calc()
end

#
function _channel_update()
end



# SISO復号器
function _siso_decoder(;Ld, decoder)
    if decoder isa TurboDecoder
        flag = false
        decoded_llrs, flag = deocder(Ld)
    elseif decoder isa LDPCDecoder
        decoded_llrs, flag = deocder(Ld)
    end
    decodedr_llrs, flag
end
