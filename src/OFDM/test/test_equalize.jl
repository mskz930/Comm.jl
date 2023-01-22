# ChannelEqualizeの動作テスト

using Random
using Commun
using Commun.Digimod, Commun.Comch, Commun.OFDM
using Commun.OFDM.Chest: convert_to_cfr

Random.seed!(1234);

n_bit = 4000; bps = 2; M = 2^bps
n_fft, n_gi, n_dims = 64, 16, (2,2)
n_data = n_bit / (M * n_dims[1])

ofdm = Ofdm(64, 16, n_dims, Pilot(:comb, df=4, t0=1))
virtualmap!(ofdm, n_data)
framelen = length(ofdm.maplist)
frame = OfdmFrame(ofdm, framelen, pilot_insert=true)

tx_bits = bitrand(1000);
tx_syms = mod(Qam(4), tx_bits)
tx_frame = copy(frame)
tx_sig  = ofdm_mod(ofdm, tx_syms, frame=tx_frame)

pdp = [1.0, 0.0, 1.0]; pdp ./= sum(pdp)
rx_sig, chresp = multipath_fading_channel(tx_sig, pdp, n_dims, fdTs=0.0)
# rx_sig = awgn_channel!(rx_sig, 0.0001)

rx_frame = ofdm_demod(ofdm, rx_sig, frame=tx_frame)
CFR = convert_to_cfr(ofdm, chresp)

demod_data = channel_equalize(ofdm, rx_frame, CFR, method=:MMSE)
