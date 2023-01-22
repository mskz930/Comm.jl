using Revise
using Commun.Comch
using Commun.Digimod
using Commun.OFDM
using Commun.OFDM.CPR: gen_coef_matrix
using Random: bitrand

# PDP
pdp = zeros(10);
pdp[:] .= 1.0;
pdp ./= sum(pdp)

qam = Qam(4)
ofdmp = Ofdm(n_fft=64, n_gi=5, n_dims=(1,1), pilot=PilotType())
ofdm_frame = OfdmFrame(ofdmp, framelen=20)

tx_frame = deepcopy(frame)

tx_bits = bitrand(1000)
tx_data = mod(qam, tx_bits)
tx_sig  = ofdm_mod(ofdmp, tx_frame, tx_data)
rx_sig, chresp = multipath_fading_channel(tx_sig, pdp, ofdmp.n_dims, fdTs=0.0)
rx_sig = awgn_channel!(rx_sig, 0.01)

Juno.@enter gen_coef_matrix(chresp, ofdmp.n_fft, ofdmp.n_gi, ofdmp.n_dims, domain=:freq)

hard_cpr_in_freq(ofdmp, rx_sig, chresp, x->slice(qam, x), n_iter=1)
