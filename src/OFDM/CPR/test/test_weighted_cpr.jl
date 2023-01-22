# test for weighted CPR in time domain

using Revise
using Commun.Comch
using Commun.Digimod
using Commun.OFDM
using Random: bitrand

# using OFDM.CPR: weighted_cpr

# パラメータ
pilotp = (type=:comb, Δf=div(64, 16), Δt=2, t0=1)
qam = Qam(4) # QPSK
ofdmp = Ofdm(n_fft=64,
             n_gi=16,
             n_dims=(1,1),
             pilot_params=pilotp
) # OFDMパラメータ

ofdm_frame = OfdmFrame(ofdmp, framelen=10, pilot_insert=true)

function trasnceiver(ofdmp, qam, tx_frame, n_bit, frame_len, pdp)
    # 送信データ
    tx_bits  = bitrand(n_bit)
    tx_syms  = mod(qam, tx_bits) # QAM変調
    tx_sig = ofdm_mod(ofdmp, tx_frame, tx_syms)

    # 通信路
    pdp = ones(5); pdp ./= length(pdp)
    rx_sig, chresp = multipath_fading_channel(tx_sig, pdp, ofdmp.n_dims, fdTs=0.0, istail=false)
    # awgn_channel!(rx_sig, 0.01)
    rx_frame = copy(tx_frame)
    rx_frame = ofdm_demod(ofdmp, rx_frame, rx_sig)

    # CFR = OFDM.ChEst.convert_to_cfr(ofdmp, chresp)
    rx_data = weighted_cpr(ofdmp, rx_frame.sig, chresp)
    # channel_equalize(ofdmp, rx_frame, CFRs)
end
