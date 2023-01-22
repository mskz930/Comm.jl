"""OFDM変調
"""
function ofdmmod(ofdm::OfdmParams, tx_frame::AbstractArray)
    # IFFT, CP付加, 並直列変換
    tx_frame   = sqrt(ofdm.n_fft) .* ifft(tx_frame, 1)
    tx_ofdm_sig = vcat(view(tx_frame, ofdm.n_fft-ofdm.n_gi+1:ofdm.n_fft,:,:), tx_frame)
    tx_ofdm_sig = reshape(ofdm_sig, :, ofdm.n_dims[1]) # 並直列変換
    tx_ofdm_sig
end

# データマッピングあり(パイロットを含む)
function ofdmmod(ofdm::OfdmParams, frame::AbstractArray, tx_data::AbstractVector)
    # シンボルマップ, IFFT, CP付加, 並直列変換
    tx_frame, frame_index = subc_mapping(ofdm, frame, tx_data)
    ifft_out = sqrt(ofdm.n_fft) .* ifft(tx_frame, 1)
    tx_ofdm_sig  = cp_add(ifft_out, ofdm.n_fft, ofdm.n_gi)
    tx_ofdm_sig  = reshape(tx_sig, :, ofdm.n_dims[1])
    tx_ofdm_sig, tx_frame, frame_index
end

# データマッピングあり(パイロットは含まない)
function ofdmmod(ofdm::OfdmParams, frame::AbstractArray, frame_index::AbstractArray, tx_data::AbstractVector)
    # シンボルマップ, IFFT, CP付加, 並直列変換
    tx_frame, frame_index = data_mapping(ofdm, frame, frame_index, tx_data)
    ifft_out = sqrt(ofdm.n_fft) .* ifft(tx_frame, 1)
    tx_sig   = cp_add(ifft_out, ofdm.n_fft, ofdm.n_gi)
    tx_sig   = reshape(tx_sig, :, ofdm.n_dims[1])
    return tx_sig, tx_frame, frame_index
end




"""
サブキャリアマッピング
"""
function subcarrier_mapping!(ofdm::OfdmParams, tx_frame, tx_symbols)
    # パラメータ
    n_syms = length(tx_symbols) # 送信シンボル数
    findex = zeros(Int64, size(tx_frame)) # frameindex

    # パイロットシンボル生成
    n_pilots = length(ofdm.indices["pilot"]) # パイロット数
    pilot_time_inds = gen_pilot_time_ids(ofdm.pilot, frame_size)
    pilot_symbols = randqam(n_pilots; M=4)

    # サブキャリアマッピング
    tx_symbols = _divide_stream(tx_symbols, ofdm.n_dims[1]) # ストリーム分割
    pilot_symbols = _split(pilot_symbols, length.(ofdm.pilot.indices)) # 配列分割
    inds = sort!(collect(union(ofdm.indices["data"], ofdm.indices["pilot"]))) # データindex
    for p in axes(tx_frame,3)
        len = length(tx_symbols[p]) # シンボル数/ストリーム
        rem = len # 未割当てのシンボル残数
        for n in axes(tx_frame,2)
            id = pilot_time_inds[n,p] # パイロットID
            flag = id > 0 ? true : false;
            pinds = flag && !isempty(ofdm.pilot.indices[1]) ? ofdm.pilot.indices[id] : Set(Int[])
            l = 1;
            for k in inds
                if !(k in pinds) && rem>0
                    tx_frame[k,n,p] = tx_symbols[p][len-rem+1]
                    findex[k,n,p] = -1 # data
                    rem -= 1
                elseif (k in pinds)
                    tx_frame[k,n,p] = pilot_symbols[id][l]
                    findex[k,n,p] = id # pilot
                    l += 1
                end
            end
            rem==0 && (frame_len=n; break)
        end
    end
    tx_frame[:,1:frame_len,:], findex[:,1:frame_len,:]
end

# サブキャリアマッピング
function subc_mapping(tx_syms, tx_frame, findex, n_dims, pilot)
    n_sym = length(tx_syms)
    tx_syms = serial_to_parallel(tx_syms, n_dims[1])
    pilot_time_inds = _get_pilot_time_inds()

    frame_data = tx_frame.data
    cnt = 0; maxj = size(frame_data,2)
    for j in axes(frame_data,2) # time
        # パイロット挿入
        for k in axes(frame_data,3) # space
            id = tx_ord[state, k]
            if pilot_insertion && id > 0
                for (i,m) in enumerate(pinds[id])
                    frame_data[i,j,k] = tx_pilots[id][m]
                    findex[i,j] = k
                end
            end
        end
        # データシンボル挿入
        for i in axes(frame_data,1)
            if findex[i,j] < 0 && cnt <= size(tx_syms,2)
                for k in axes(frame_data,3)
                    frame_data[i,j,k] = tx_syms[k, cnt]
                    cnt += 1
                end
            end
        end
        # 終了条件
        (cns > size(tx_data,2)) && break
        maxj = j
    end
    # フレームに情報を書き込む
    tx_frame.len = maxj
    tx_frame.n_sym = n_sym

    tx_data, tx_frame
end


# データシンボルをマッピングする
function data_mapping!(frame::AbstractArray{<:Number}, tx_symbols::AbstractVector, frameindex)
    len = length(tx_symbols) # シンボル数
    n_tx = size(frame)[end] # フレームの外側次元
    tx_symbols = _split(tx_symbols, n_tx) # ストリーム分割
    for p in axes(frame,3) # 空間軸
        for n in axes(frame,2)
            len = length(tx_symbols[p])
            rem = len
            for k in axes(frame,1)
                if frameindex[i,j,k] < 0
                    if rem > 0
                        frame[i,j,k] = tx_symbols[p][len-rem+1]
                        rem -= 1
                    end
                end
            end
            rem == 0 && break
        end
    end
end




"""
    data_demapping!(ofdm::OfdmParams, rx_sig)

データシンボルの抽出
"""
function data_demapping(ofdm::OfdmParams, rx_ofdmsig)
    N = ofdm.n_fft # FFT数
    rx_ofdmsig = reshape(rx_ofdmsig, ofdm.n_fft+ofdm.n_gi, :, ofdm.n_dims[2]) # 直並列変換
    rx_frame = rx_ofdmsig[ofdm.n_gi+1:size(rx_ofdmsig,1), :, :] # CP除去
    rx_frame = fft!(rx_frame, 1) ./ sqrt(N) # FFT
    # rx_frame = rx_frame[order,:,:] # 周波数の正負入れ替え
end
