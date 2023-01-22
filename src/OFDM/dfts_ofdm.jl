"""
SC-FDM変調
"""
function scfdm_mod(ofdm::Ofdm, tx_frame::AbstractArray{T,3}, tx_data) where T <: Complex
    @unpack n_fft, n_gi, n_dim = ofdm
    # DFT拡散マッピング
    tx_frame = dfts_mapping!(params, pilot, copy(data_frame), tx_symbol, frameinfo)

    # IFFT
    tx_frame = sqrt(N) .* ifft(tx_frame, 1)
    tx_frame = gi_add(tx_frame, n_fft, n_gi)# CP付加
    tx_sig = reshape(tx_frame, :, n_dim[1])
    tx_sig
end

# DFT拡散マッピング
function dft_s_mapping(ofdm, frame, tx_data; M=length(ofdm.idxs[:data]))
    tx_frame = copy(frame)
    didxs = ofdm.idxs[:data]
    tx_data = serial_to_parallel(tx_data, M, :zero)
    tx_data = fft(tx_data, 1) ./ sqrt(M)
    tx_data = vec(tx_data)
    tx_frame = subcmap(ofdm, frame, tx_data)
    tx_frame
end

function dft_s_demapping(ofdm, rx_data, M=length(ofdm.idxs[:data]))
    rx_data = serial_to_parallel(rx_data, M, :zero)
    rx_data = ifft!(rx_data, 1) .* sqrt(M)
    return rx_data
end
