# OFDM用自作ツール関数


"""送信フレーム配列から有効な送信データシンボルを取り出す
"""
function extract_data(ofdm::Ofdm, frame; dims=1)
    didxs = get_idxs(ofdm.frame, :data)
    ndata = ofdm.frame.on[:data] # 有効な送信データシンボル数
    if dims==1
        dataseq = frame[:,didxs]
    else
        dataseq = frame[didxs, :]
    end
    return @view dataseq[1:ndata]
end

""" フレームからパイロットシンボルを取り出す
"""
function extract_pilot(ofdm::Ofdm, frame)
    tidxs = get_idxs(ofdm.frame, :pilot) # timeidxs
    return @view frame[pidxs, :]
end

# ガードキャリアから雑音電力を推定する
function nvar_estimate(ofdm, rx_frame)
    gidxs = ofdm.idxs[:guard]
    isempty(gidxs) && error("ガードキャリアが割り当てられていません．")
    n_fft, nsym, nrx = size(rx_frame)
    guards = @view rx_frame[gidxs, :, :]
    nvar = mean(abs2, guards, dims=[1,2]) # 時間周波数平均
    return vec(nvar)
end
