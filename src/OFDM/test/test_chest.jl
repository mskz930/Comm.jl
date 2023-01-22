using Random
using Commun.Digimod
using Commun.Comch
using Commun.OFDM

pdp = [1.0, 0.0, 0.0, 0.0, 1.0]
pdp ./= sum(pdp)

SNR = 30
N0  = 10^(-SNR/10)

modtype = Qam(4)

n_fft  = 64
n_gi   = 16
n_dims = (1,1)


ofdm  = Ofdm(n_fft, n_gi, n_dims, Pilot(:comb, t0=1, dt=2, df=8))
frame = OfdmFrame(ofdm, framelen=20, pilot_insert=true)

tx_bits = bitrand(10000)
tx_data = mod(modtype, tx_bits)
tx_frame = deepcopy(frame)
tx_sig = ofdm_mod(ofdm, tx_frame, tx_data)


rx_sig, chresp = multipath_fading_channel(tx_sig, pdp, n_dims)
rx_sig = awgn_channel!(rx_sig, N0)

rx_frame = ofdm_demod(ofdm, tx_frame, rx_sig)


CFR = channel_estimate(ofdm, rx_frame, tx_frame, average=true)
