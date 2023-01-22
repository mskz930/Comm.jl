include("../../Ofdm.jl")
include("../../../../Channel/Channel.jl")
nbits = 2000
M = 4
qam = Ofdm.QamMod(M)
nfft=64; cpsize=16; npi=12; ngc=6; ndim=(2,2)
params = Ofdm.OfdmParams(nfft,cpsize,npi,ngc,ndim; pilot_type=:comb, pilot_space=4, pilot_interval=2)
indices = Ofdm.get_indices(params)
frame = Ofdm.gen_frame(params, indices, 20)
pdp = Channel.exponent(10,20)
# 送信機
tx_bit = rand(nbits) .> 0.5
tx_symbol = Ofdm.qammod(qam, tx_bit)
tx_ofdmsig = Ofdm.ofdmmod(params, frame, tx_symbol, indices)

# 通信路
rx_ofdmsig = Channel.multipath_fading(tx_ofdmsig, ndim, pdp)
