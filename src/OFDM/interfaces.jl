# interfaces.jl

# 送信機
function transmitter(ofdm::Ofdm, nbit, mod_t)
    tx_bits = randbit(nbit)
    tx_data = mod(modtype, tx_bits)
    tx_frame = subcmap(ofdm, frame, tx_data)
    tx_sig = mod(ofmd, tx_frame)
    tx_sig
end


# 受信機
function recevier(ofdm::Ofdm, mod_t, nbit, frame, N0, chresp)
    rx_frame = demod(ofdm, rx_sig)
    cfr = to_cfr(ofdm, chresp)
    rx_data  = equalize(ofdm, rx_frame, cfr)
    rx_bits = demod(mod_t, rx_data)
    return rx_bits
end


